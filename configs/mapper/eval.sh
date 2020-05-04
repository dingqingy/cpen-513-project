cd ../../
scons -j4 --accelergy --d
cd -
../../build/timeloop-mapper VGGN.yaml ERT.yaml > out.log
# gdb --args ../../build/timeloop-mapper VGGN.yaml ERT.yaml

# echo 'clustered balance'
# ../../build/timeloop-mapper LB_CONV1_CK.yaml > CK_CLB.log
# 
# echo 'global balance'
# export TIMELOOP_GLOBAL_SORT=True
# ../../build/timeloop-mapper LB_CONV1_CK.yaml > CK_GLB.log
# unset TIMELOOP_GLOBAL_SORT
# echo 'sync PE'
# export TIMELOOP_USE_SYNC_PE=True
# ../../build/timeloop-mapper LB_CONV1_CK.yaml > CK_NLB.log
# unset TIMELOOP_USE_SYNC_PE
# ../../build/timeloop-mapper PQ.yaml > PQ.log
# ../../build/timeloop-mapper full-N.yaml > full-N.log
# ../../build/timeloop-mapper FC2.yaml > FC2.log
# ../../build/timeloop-mapper densePQ.yaml > densePQ.log
# gdb --args ../../build/timeloop-mapper FC2.yaml

# CK partition
# ../../build/timeloop-mapper VGG.cfg > VGG-32bit.log
# gdb --args ../../build/timeloop-mapper VGG.cfg

# CN partition
# ../../build/timeloop-mapper VGGN.cfg > VGGN-32bit.log
# ../../build/timeloop-mapper VGGN_bk.cfg > VGGN-32bit_bk.log
# gdb --args ../../build/timeloop-mapper VGGN.cfg
# cp timeloop-mapper.map.cfg ../evaluator/superblock/bwK8.cfg
