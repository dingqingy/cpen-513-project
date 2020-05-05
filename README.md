# CPEN 513 Project
Optimizing CNNAccelerator Dataflow Using Simulated Annealing

## Report
* [Report](project_report.pdf)

## Getting started
Follow instruction in timeloop_setup.md, make sure to have all dependencies and other setup
Run simulated annealing mapping search algorithm on an eyeriss like hardware on a VGG_conv_5_2 layer
```
cd configs/mapper/
bash eval.sh
```

Run tuning on more benchmarks layers from resnet18
```
cd scripts/
python sample.py
```
Note that the current script is not capable of handling the error caused from the case that the primitive list is initialized to a bad position where all mapping produced is rejected, and thus the search fails with no output generated. The result parsing script will look for the output and if it doesn't find it, it produces the error. This error handling will be provided in the future (or we change the code such that reinitialize if landed in a bad region in the mapspace)