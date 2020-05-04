#pragma once

#include <iterator>
#include <mutex>
#include <regex>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>

#include "util/numeric.hpp"
#include "util/misc.hpp"
#include "workload/problem-shape.hpp"
#include "mapspaces/mapspace-base.hpp"
#include "mapspaces/subspaces.hpp"
#include "compound-config/compound-config.hpp"
#include "mapping/arch-properties.hpp"
#include "mapping/constraints.hpp"

namespace mapspace
{

//--------------------------------------------//
//                Didi MapSpace               //
//--------------------------------------------//

// struct Primitive
// {
// 	problem::Shape::DimensionID id;
//   std::string name;
// 	unsigned val;

// 	void Print() const
// 	{
// 		std::cout <<name<<val<<" ";
// 	}
// };

typedef std::pair<problem::Shape::DimensionID,unsigned> Primitive;


class Didi : public MapSpace
{
 protected:
  uint128_t timestamp_;
 	std::vector<Primitive> primitive_list_;
  std::vector<Primitive> proposed_; // we evaluate on proposed, update if we take the swap

 	// split is disabled given the scope of the project
 	std::vector<Didi*> splits_; // always have size 1

 	// Abstract representation of the architecture.
  ArchProperties arch_props_;

  // Constraints.
  mapping::Constraints constraints_;

 public:
  Didi(
    config::CompoundConfigNode config,
    config::CompoundConfigNode arch_constraints,
    model::Engine::Specs arch_specs,
    const problem::Workload& workload,
    bool skip_init = false) :
      MapSpace(arch_specs, workload),
      arch_props_(arch_specs),
      constraints_(arch_props_, workload)
  {
    if (!skip_init)
    {
      Init(config, arch_constraints);
    }
  }

  //------------------------------------------//
  //        Initialization and Setup          // 
  //------------------------------------------//
  
  //
  // Init() - called by derived classes or by constructor.
  //
  void Init(config::CompoundConfigNode config, config::CompoundConfigNode arch_constraints)
  {
    // Parse user config.
    Parse(config, arch_constraints);

    InitPrimitivePairs();
    timestamp_=0;

    // Sanity check
    std::cout<< "Checking Init Primitive list"<<std::endl;
    for(unsigned i = 0; i < primitive_list_.size(); i++)
    {
    	// primitive_list_[i].Print()
      std::cout << problem::GetShape()->DimensionIDToName.at(primitive_list_[i].first) << primitive_list_[i].second << " " << std::endl;
    }
    std::cout << std::endl;

    std::cout<< "Checking Init Proposed Primitive list"<<std::endl;
    for(unsigned i = 0; i < proposed_.size(); i++)
    {
      // proposed_[i].Print()
      std::cout << problem::GetShape()->DimensionIDToName.at(proposed_[i].first) << proposed_[i].second << " " << std::endl;
    }
    std::cout << std::endl;
  }

  void InitPrimitivePairs()
  {
  	// Note: User factor should not be considered for swap only within the tiling level
  	// to simplify: assume no user factor (discard) first?
    // TODO: fix it to allow user factor (maybe first allow spatial and then arbitrary constraints)
  	auto user_factors = constraints_.Factors();
    auto user_max_factors = constraints_.MaxFactors();

    assert(user_factors.size() <= arch_props_.TilingLevels());

    // Done: based on the workload, produce primitive_list in an arbitrary order
    // TODO: we can fix spatial level to reduce the complexity

    primitive_list_.clear();
    for (unsigned idim = 0; idim < unsigned(problem::GetShape()->NumDimensions); idim++) // idim for dim id
    {
      auto bound = workload_.GetBound(idim);
      while (bound % 2 == 0)  
      {  
        primitive_list_.push_back(std::make_pair(idim, 2));
        bound /= 2;  
      }

      for (int i = 3; i <= std::sqrt(bound); i += 2)  
      {  
        while (bound % i == 0)  
        {  
          primitive_list_.push_back(std::make_pair(idim, i));
          bound = bound/i;  
        }  
      }

      if(bound>2)
      {
        primitive_list_.push_back(std::make_pair(idim, bound));
      }
    }

    // Done: add special premitives such as level seperators    
    for (unsigned i = 0; i<arch_props_.TilingLevels()-1;i++)
    {
      primitive_list_.push_back(std::make_pair(problem::GetShape()->NumDimensions, 0)); // level seperator
    }    

    std::shuffle(std::begin(primitive_list_), std::end(primitive_list_), std::default_random_engine());

    // TODO: infer spatial primitives implicitly based on level seperator and IsSpatial(level)

    proposed_ = primitive_list_;
  }

  Didi(const Didi& other) = default;

  // Done: Split should produce no effect in simulated anneadling 
  // Or in the mapper, disable split if use SA search

  std::vector<MapSpace*> Split(std::uint64_t num_splits)
  {
    num_splits = 1;

    std::cout << "No parallelization with SA, num_splits set to "<< num_splits << std::endl;

    std::vector<Didi*> splits;
    std::vector<MapSpace*> retval;
    for (unsigned i = 0; i < num_splits; i++)
    {
      Didi* mapspace = new Didi(*this);
      
      splits.push_back(mapspace);
      retval.push_back(static_cast<MapSpace*>(mapspace));
    }

    splits_ = splits;
    return retval;
  }

  // Done: InitPruned should be disabled
  void InitPruned(uint128_t index_factorization_id)
  {
  	(void) index_factorization_id;
  }

  //------------------------------------------//
  //           Mapping Construction           // 
  //------------------------------------------//
  
  //
  // ConstructMapping()
  //   Given a primitive list representation this map space,
  //   construct a full mapping.
  //

  bool ConstructMapping(
      mapspace::ID mapping_id, // what I really want is timestamp
      Mapping* mapping)
  {
  	(void) mapping_id;

  	// A set of subnests, one for each tiling level.
    loop::NestConfig subnests(arch_props_.TilingLevels());

    // Note in all stage subnests modified inplace
    // // === Stage 0 ===
    // InitSubnests(subnests); 

    // // === Stage 1 ===
    // PermuteSubnests(subnests); // reorder based on proposed primitive list

    // // === Stage 2 ===
    // AssignIndexFactors(subnests); // merge primitive with the same name
    PrimitiveToSubNest(subnests);

    // === Stage 4 ===
    mapping->datatype_bypass_nest = ConstructDatatypeBypassNest(); // simplify this with assuming no bypass

    // FIXME: optimization: if we are using the mapspace in deferred/pruned mode then
    // the index factorization does not need to be re-processed.

    // We had to reverse the order of stage 4 and 3 because AssignSpatialTilingDirections
    // needs the datatype bypass nest to determine if a spatial fanout is possible or
    // not.
    
    // === Stage 3 ===
    bool success = AssignSpatialTilingDirections(subnests, mapping->datatype_bypass_nest);
    if (!success)
    {
      return false;
    }

    // Check possible failure reason for Uber AssignSpatialTilingDirections
    // If we choose to fix spatial here ???? (think again)

    // Concatenate the subnests to form the final mapping nest.    
    std::uint64_t storage_level = 0;
    for (uint64_t i = 0; i < arch_props_.TilingLevels(); i++)
    {
      uint64_t num_subnests_added = 0;
      for (unsigned dim = 0; dim < subnests[i].size(); dim++) // Note: already pruned on each tile level
      {
        // Ignore trivial factors
        // This reduces computation time by 1.5x on average.
        if (subnests[i][dim].start + subnests[i][dim].stride < subnests[i][dim].end)
        {
          mapping->loop_nest.AddLoop(subnests[i][dim]);
          num_subnests_added++;
        }
      }
      if (!arch_props_.IsSpatial(i))
      {
        if (num_subnests_added == 0)
        {
          // Add a trivial temporal nest to make sure
          // we have at least one subnest in each level.
          mapping->loop_nest.AddLoop(problem::Shape::DimensionID(int(problem::GetShape()->NumDimensions) - 1),
                                     0, 1, 1, spacetime::Dimension::Time);
        }
        mapping->loop_nest.AddStorageTilingBoundary();
        storage_level++;
      }
    }

    // Finalize mapping.
    // mapping->id = mapping_id.Integer();
    // Done: figure out timestamp increase, should happen in search or here?
    // we do it in search, here just use another timestamp_ to distinguish mapping
    mapping->id = timestamp_;
    timestamp_++;
    return true;
  }

  //
  // Mapping Construction
  // corresponds to Uber Stage 0-2
  //
  void PrimitiveToSubNest(loop::NestConfig& subnests)
  {
    unsigned start = 0;
    unsigned end = 0;
    std::vector<unsigned> merged_bounds;
    std::vector<unsigned> primitive_counts;
    std::vector<double> primitive_relative_pos; // sum
    // std::vector<std::pair<unsigned, unsigned>> order_arbitrator; // track primitive appearance count and sum of the relative position to start
    // for (int)

    for (uint64_t level = 0; level < arch_props_.TilingLevels(); level++)
    {
      // identify position of level seperator
      while(end < proposed_.size() || IsLevelSeperator(proposed_[end]))
      {
        end++;
      }

      // init auxiliary structure
      merged_bounds.clear();
      merged_bounds.resize(problem::GetShape()->NumDimensions, 1);

      primitive_counts.clear();
      primitive_counts.resize(problem::GetShape()->NumDimensions, 0);

      primitive_relative_pos.clear();
      primitive_relative_pos.resize(problem::GetShape()->NumDimensions, 0);

      // order_arbitrator.clear();
      // for (i = 0; i < problem::GetShape()->NumDimensions;i++)
      //   order_arbitrator.push_back(std::make_pair(0, 0));

      // merge primitives and determine the order
      for(unsigned i = start; i<end; ++i)
      {
        dim = proposed_.at(i).first;
        val = proposed_.at(i).second;

        merged_bounds.at(dim) *= val;
        primitive_counts.at(dim) += 1;
        primitive_relative_pos += i - start;
        // order_arbitrator.at(dim).first++;
        // order_arbitrator.at(dim).second += i-start;
      }

      std::vector<std::pair<double, problem::Shape::DimensionID>> arbitrator; // arbitrator has relative order and dimID
      arbitrator.clear();
      for (int idim = 0; idim < int(problem::GetShape()->NumDimensions); idim++)
      {
        if(primitive_counts[idim]>0)
        {
          arbitrator.push_back(std::make_pair(primitive_relative_pos[idim]/primitive_counts[idim], idim))
        }
      }

      if(!arbitrator.empty())
        std::sort(arbitrator.begin(), arbitrator.end());

      for(unsigned i = 0; i < arbitrator.size(); i++)
      {
        auto spacetime_dim = arch_props_.IsSpatial(level)
          ? spacetime::Dimension::SpaceX // Placeholder.
          : spacetime::Dimension::Time;
        
        
        auto dim = arbitrator[i].second;
        loop::Descriptor loop;
        loop.dimension = problem::Shape::DimensionID(dim);
        loop.start = 0;
        loop.end = merged_bounds[dim];
        loop.stride = 1;                           // FIXME.
        loop.spacetime_dimension = spacetime_dim;
        
        subnests.at(level).push_back(loop);
      }

      start = end + 1;
    }
    // TODO: missing handling for spatial level and no allocated temporal dim
    // TODO: check how is tiling boundary expressed
  }

  bool IsLevelSeperator(Primitive p) const
  {
    return (p.first == problem::GetShape()->NumDimensions) && (p.second == 0)
  }

  //
  // Mapping Construction
  // Stage 4: Construct datatype bypass nest. we assume no bypass for simplicity
  //
  tiling::CompoundMaskNest ConstructDatatypeBypassNest()
  {
    tiling::CompoundMaskNest seed_mask_nest;
    for (unsigned pvi = 0; pvi < unsigned(problem::GetShape()->NumDataSpaces); pvi++)
    {
      for (unsigned level = 0; level < arch_specs_.topology.NumStorageLevels(); level++)
      {
        seed_mask_nest.at(pvi).set(level);
      }
    }
    return seed_mask_nest;
  }
/*
  //
  // Mapping Construction
  // Stage 0: Initialize subnests.
  //
  void InitSubnests(loop::NestConfig& subnests)
  {
    // Construct num_storage_levels loop-nest partitions and assign dimensions.
    // This is the only stage at which the invariant subnests[][dim].dimension == dim
    // will hold. The subnests will later get permuted, breaking the invariant.
    for (uint64_t level = 0; level < arch_props_.TilingLevels(); level++)
    {
      auto spacetime_dim = arch_props_.IsSpatial(level)
        ? spacetime::Dimension::SpaceX // Placeholder.
        : spacetime::Dimension::Time;
        
      // Each partition has problem::GetShape()->NumDimensions loops.
      for (int idim = 0; idim < int(problem::GetShape()->NumDimensions); idim++)
      {
        loop::Descriptor loop;
        loop.dimension = problem::Shape::DimensionID(idim); // Placeholder.
        loop.start = 0;
        loop.end = 0;                              // Placeholder.
        loop.stride = 1;                           // FIXME.
        loop.spacetime_dimension = spacetime_dim;
        
        subnests.at(level).push_back(loop);
      }
    }
  }

  //
  // Mapping Construction
  // Stage 1: Permute Subnests.
  //
  void PermuteSubnests(loop::NestConfig& subnests) // another argument to indicate how to swap
  {
    loop::NestConfig reordered(arch_props_.TilingLevels());
    
    // Obtain a pattern of loop variables for all levels.
    // auto dimensions = permutation_space_.GetPatterns(mapping_permutation_id);
    // assert(dimensions.size() == subnests.size());

    for (uint64_t level = 0; level < arch_props_.TilingLevels(); level++)
    {
      // Re-order the subnest based on the pattern. 
      // assert(dimensions[level].size() == subnests[level].size());
      // for (unsigned i = 0; i < dimensions[level].size(); i++)
      // {
      //   auto target_dim = dimensions[level][i];
      //   assert(subnests[level][int(target_dim)].dimension == target_dim);
      //   reordered[level].push_back(subnests[level][int(target_dim)]);
      // }

      // TODO: perform the swap and legalize the loop order (deduplication)
            // Each partition has problem::GetShape()->NumDimensions loops.
      for (int idim = 0; idim < int(problem::GetShape()->NumDimensions); idim++)
      {
        loop::Descriptor loop;
        loop.dimension = problem::Shape::DimensionID(idim); // Placeholder.
        loop.start = 0;
        loop.end = 0;                              // Placeholder.
        loop.stride = 1;                           // FIXME.
        loop.spacetime_dimension = spacetime_dim;
        
        subnests.at(level).push_back(loop);
      }
    }

    subnests = reordered;
  }

  //
  // Mapping Construction
  // Stage 2: Assign Index Factors.
  //
  void AssignIndexFactors(loop::NestConfig& subnests) // another argument to indicate how to swap? no need if already legalized
  {
    for (uint64_t level = 0; level < arch_props_.TilingLevels(); level++)
    {
      // for (auto& loop : subnests[level])
      // {
      //   loop.end = int(index_factorization_space_.GetFactor(
      //                    mapping_index_factorization_id,
      //                    loop.dimension,
      //                    level));
      // }

      // TODO: merge assigned index
    }
  }
*/
  // validate spatial fanout
  // TODO: get spatial primitive to work properly

  //
  // Mapping Construction
  // Stage 3: Decide which of the spatial loop nests are along the space_x dimension.
  //
  bool AssignSpatialTilingDirections(loop::NestConfig& subnests,
                                     tiling::CompoundMaskNest datatype_bypass_nest)
  {
    (void) datatype_bypass_nest;
    bool success = true;

    // TODO: implicitly infer split
    // auto spatial_splits = spatial_split_space_.GetSplits(mapping_spatial_id);
    //auto datatype_bypass_masks = tiling::TransposeMasks(datatype_bypass_nest);
    
    double cumulative_fanout_utilization = 1.0;

    for (uint64_t level = 0; level < arch_props_.TilingLevels() && success; level++)
    {
      if (!arch_props_.IsSpatial(level))
      {
        continue;
      }
      
      // Note that spatial levels are never bypassed. Therefore, we do not have
      // to deal with the bypass mask.
      // auto& datatype_bypass_mask = datatype_bypass_masks.at(storage_level-1);

      success &= AssignSpatialTilingDirections_Level_Expand(
        // spatial_splits.at(level),
        subnests[level],
        level,
        cumulative_fanout_utilization);
      
    } // for (level)
    
    // success &= (cumulative_fanout_utilization >= constraints_.MinParallelism()); // ignored for simplicity
      
    return success;
  }

  bool AssignSpatialTilingDirections_Level_Expand(// std::uint32_t spatial_split,
                                                  std::vector<loop::Descriptor>& level_nest,
                                                  unsigned tiling_level_id,
                                                  double& fanout_utilization)
  {
    // This version of the function assumes that spatial tiling will expand
    // the instances for *each* datatype exactly by the tiling parameters. For
    // example, if K=16 is a spatial factor, then 16 instances of the next
    // inner level will be created for Weights, Inputs and Outputs.
    
    bool success = true;

    unsigned storage_level_id = arch_props_.TilingToStorage(tiling_level_id);
    auto level_specs = arch_specs_.topology.GetStorageLevel(storage_level_id);

    std::size_t x_expansion = 1;
    std::size_t y_expansion = 1;
    
    // Based on the spatial mapping ID, split the level nest into two sections:
    // first X and then Y.
    for (unsigned i = 0; i < level_nest.size(); i++)
    {
      auto& loop = level_nest.at(i);
      
      assert(loop::IsSpatial(loop.spacetime_dimension));
      assert(loop.stride == 1);

      if (x_expansion * (loop.end - loop.start) <= arch_props_.FanoutX(storage_level_id)) // is x possible to expand
      {
        // X
        x_expansion *= (loop.end - loop.start);
        loop.spacetime_dimension = spacetime::Dimension::SpaceX;
      }
      else if (y_expansion * (loop.end - loop.start) <= arch_props_.FanoutY(storage_level_id)) // is y possible to expand
      {
        // Y
        y_expansion *= (loop.end - loop.start);
        loop.spacetime_dimension = spacetime::Dimension::SpaceY;
      }
      else
      {
        success = false;
      }
    }

    std::size_t fanout_max;
    
    // if (level_specs->SharingType() == model::DataSpaceIDSharing::Shared)
    // {
    // if (x_expansion > arch_props_.FanoutX(storage_level_id))
    //   success = false;
      
    // if (y_expansion > arch_props_.FanoutY(storage_level_id))
    //   success = false;

    // fanout_max = arch_props_.Fanout(storage_level_id);
    // }
    // else
    // {
    //   std::size_t x_fanout_max = 0;
    //   std::size_t y_fanout_max = 0;

    //   // The following loop is silly since we now only allow one fanout per level
    //   // (as opposed to a per-dataspace fanout for partitioned levels). However,
    //   // we will keep the code because we may need to move to a multiple-buffers
    //   // per level later. The loop will not be over data spaces but buffer
    //   // instances per level.

    //   for (unsigned pvi = 0; pvi < unsigned(problem::GetShape()->NumDataSpaces); pvi++)
    //   {
    //     // auto pv = problem::Shape::DataSpaceID(pvi);

    //     if (x_expansion > arch_props_.FanoutX(storage_level_id))
    //       success = false;

    //     if (y_expansion > arch_props_.FanoutY(storage_level_id))
    //       success = false;

    //     // Track max available (not utilized) fanout across all datatypes.
    //     x_fanout_max = std::max(x_fanout_max, arch_props_.FanoutX(storage_level_id));
    //     y_fanout_max = std::max(y_fanout_max, arch_props_.FanoutY(storage_level_id));
    //   }

    //   fanout_max = x_fanout_max * y_fanout_max;
    // }

    // Compute fanout utilization at this level.
    // Ignore bypass and partitioning. The only purpose of this is to accumulate
    // the level-wise utilizations to compute arithmetic utilization.
    // fanout_utilization *= double(x_expansion) * double(y_expansion) / fanout_max;

    // if (!success)
    // {
    //   std::cerr << "Level: " << arch_props_.StorageLevelName(storage_level_id) << std::endl;
    //   std::cerr << "  X: ";
    //   std::cerr << " expansion = " << x_expansion << " fanout = " << arch_props_.FanoutX(storage_level_id) << std::endl;
    //   std::cerr << "  Y: ";
    //   std::cerr << " expansion = " << y_expansion << " fanout = " << arch_props_.FanoutY(storage_level_id) << std::endl;
    //   std::cerr << "  util = " << fanout_utilization << std::endl;
    //   std::cerr << std::endl;
    // }
    
    return success;
  }
 

  //------------------------------------------//
  //   Manipulating Primitive representation  // 
  //------------------------------------------//

  void ProposeToSwap(unsigned p1, unsigned p2)
  {
    // proposal is based on the true representation
    proposed_ = primitive_list_;
    std::swap(proposed_.at(p1), proposed_.at(p2));
  }

  void AcceptProposal()
  {
    // note if proposal not accepted, will not update primitive list
    assert(primitive_list_.size() == proposed_.size());
    primitive_list_ = proposed_;
  }

  std::size_t GetPrimitiveCount()
  {
    // use this to make sure proposed position within vector size
    assert(primitive_list_.size() == proposed_.size());
    return primitive_list_.size();
  }


  //------------------------------------------//
  //                 Parsing                  // 
  //------------------------------------------//
  
  //
  // Parse.
  //
  void Parse(config::CompoundConfigNode config, config::CompoundConfigNode arch_constraints)
  {
    // Parse constraints.
    // We accept mapspace config and arch_constraints as separate configuration
    // trees, but as far as parsing is concerned we handle them in exactly the
    // same way. The underlying parsing methods are built to handle conflicts.
    constraints_.Parse(config);
    constraints_.Parse(arch_constraints);
  }
};


} // namespace mapspace
