# tree2sexaquark

ROOT-based standalone application for the Sexaquark search in ALICE.

## Requirements

* VC
* ROOT
* KFParticle

## Building

```
mkdir build && cd build
cmake ../ -DCMAKE_EXPORT_COMPILE_COMMANDS=1
```

## Debugging

```
valgrind --leak-check=full --suppressions=$(root-config --etcdir)/valgrind-root.supp APP &> LOG
```
