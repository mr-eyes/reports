# Indexing trial using the pure C++ kProcessor version


## Cloning & building

```bash
cd cpp_exp/
git clone git@github.com:dib-lab/kProcessor.git
cd kProcessor
git checkout batchQuery
git submodule update --init --recursive
cd ..
mkdir build && cd build
cmake ..
make
```

## Run

```bash
cd cpp_exp/build
./try_idx reps_unitigs_SRR11015356_k75.fa
```

## Result

**The huge number of colors caused slowness in insersion and inflation in memory**

<details>
  <summary>**Click to the expand the debugging log**</summary>
  
```bash
[DEBUG:: ] chunk: 1
tagsMapsize: 2528241
---------------------------
[DEBUG:: ] chunk: 2
tagsMapsize: 2608262
---------------------------
[DEBUG:: ] chunk: 3
tagsMapsize: 2777606
---------------------------
[DEBUG:: ] chunk: 4
tagsMapsize: 2977908
---------------------------
[DEBUG:: ] chunk: 5
tagsMapsize: 3161182
---------------------------
[DEBUG:: ] chunk: 6
tagsMapsize: 3414790
---------------------------
[DEBUG:: ] chunk: 7
tagsMapsize: 3720254
---------------------------
[DEBUG:: ] chunk: 8
tagsMapsize: 4228800
---------------------------
[DEBUG:: ] chunk: 9
tagsMapsize: 4793336
---------------------------
[DEBUG:: ] chunk: 10
tagsMapsize: 5384756
---------------------------
[DEBUG:: ] chunk: 11
tagsMapsize: 6016696
---------------------------
[DEBUG:: ] chunk: 12
tagsMapsize: 6680597
---------------------------
[DEBUG:: ] chunk: 13
tagsMapsize: 7370437
---------------------------
[DEBUG:: ] chunk: 14
tagsMapsize: 8095047
---------------------------
[DEBUG:: ] chunk: 15
tagsMapsize: 8846856
```
</details>.