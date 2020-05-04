#pragma once

#include <iterator>
#include <unordered_set>
#include <cmath>
#include <random>

#include "mapping/mapping.hpp"
#include "mapspaces/mapspace-base.hpp"
#include "util/misc.hpp"
#include "search/search.hpp"

namespace search
{

class SimulatedAnnealing : public SearchAlgorithm
{
 private:
  enum class State
  {
    Ready,
    WaitingForStatus,
    Terminated
  };
  
 private:
  // Config.
  mapspace::MapSpace* mapspace_;
  // Live state.
  State state_;
  uint128_t valid_mappings_;
  std::uint64_t eval_fail_count_;

  // simulated annealing states
  uint128_t timestamp_;
  double temp_;
  double previous_cost_; // used to calculate delta
  double best_cost_; 

  uint128_t max_iter_;
  unsigned early_stop_iter_;
  unsigned cooling_iter_;
  double beta_;

  // random nunber stuff
  std::random_device rd_;  //Will be used to obtain a seed for the random number engine
  std::mt19937 rng_; //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis_;

  // Done: permutation list for proposing swap
  std::vector<unsigned> perm_;

 public:
  SimulatedAnnealing(config::CompoundConfigNode config, mapspace::MapSpace* mapspace) :
      SearchAlgorithm(),
      mapspace_(mapspace),
      state_(State::Ready),
      valid_mappings_(0),
      eval_fail_count_(0),
      timestamp_(0),
      previous_cost_(std::numeric_limits<double>::max()),
      best_cost_(std::numeric_limits<double>::max()),
      beta_(0.9),
      rng_(rd_()),
      dis_(std::uniform_real_distribution<>(0.0, 1.0))
  {
    (void) config;
    // TODO: get temp_, max_iter_, early_stop_iter_, cooling_iter_, beta_ (optionally)
    // from config
    temp_ = 1000;
    max_iter_ = 1000000;
    early_stop_iter_ = 10000;
    cooling_iter_ = 10;
    beta_ = 0.9;

    //
    perm_.clear();
    for(int i = 0; i <mapspace_->GetPrimitiveCount();i++)
      perm_.push_back(i);
  }

  // Done: figure out when to propose swap, when to actually swap
  // need to review mapper
  // we manipulate in proposed list in Next
  // take the actual swap in Report

  bool Next(mapspace::ID& mapping_id)
  {
    (void) mapping_id;

    // need to propose a swap
    // and how to get the swap?
    // modify the primitive representation inplace here?

    if (state_ == State::Terminated)
    {
      return false;
    }

    assert(state_ == State::Ready);

    // mapping_id = mapspace::ID(mapspace_->AllSizes());
    // for (unsigned i = 0; i < unsigned(mapspace::Dimension::Num); i++)
    // {
    //   mapping_id.Set(i, iterator_[i]);
    // }

    // Done: propose a new swap (not through new mapping ID but through update proposed_ in the Didi mapspace)
    std::shuffle(std::begin(perm_), std::end(perm_), rng_);
    mapspace_->ProposeToSwap(perm_.at(0), perm.at(1));

    // increase the timestamp similar to update in mapping_id
    timestamp_++;


    state_ = State::WaitingForStatus;
    
    return true;
  }

  void Report(Status status, double cost = 0)
  {
    (void) cost;
    
    assert(state_ == State::WaitingForStatus);

    // bool skip_datatype_bypass = false;
    if (status == Status::Success)
    {
      valid_mappings_++;

      auto r = dis_(rng_);
      auto d_cost = cost - previous_cost_;

      if (r < std::exp(-d_cost / temp_))
      {
        // Done: swap that actually change the primitive_list from proposed_ in the Didi mapspace
        mapspace_->AcceptProposal();
        previous_cost_ = cost;
      }

      best_cost_ = std::min(best_cost_, cost);
    }
    else if (status == Status::MappingConstructionFailure)
    {
      // Accelerate search by invalidating bad spaces.
      // ConstructMapping failure =>
      //   Combination of (IF, LP, S) is bad.
      //   Skip all DBs.
      // skip_datatype_bypass = true;
    }
    else if (status == Status::EvalFailure)
    {
      // PreEval/Eval failure (capacity) =>
      //   Combination of (IF, DB) is bad.
      //   If all DBs cause Eval failure for an IF, then that IF is bad,
      //   no need to look at other LP, S combinations.
      eval_fail_count_++;
    }

    bool mapspace_remaining = timestamp_ < max_iter_; // keep searching until the termination condition
    // TODO: add early termination

    if (mapspace_remaining) //  && valid_mappings_ < search_size_)
    {
       state_ = State::Ready;
    }
    else
    {
      state_ = State::Terminated;
    }
	}
};

} // namespace search
