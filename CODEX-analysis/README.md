# CODEX data analysis
Codes for CODEX data analysis:

## Denpendency
- [deepcell](https://github.com/vanvalenlab/deepcell-tf)
- [napari](https://napari.org/) (For visualization)

## Cell segmentation
Set the parameters in `1_segment.py` and run the following command:
```bash
python 1_segment.py
```

## Quantification to get cell by gene matrix

Set the parameters in `2_quantify.py` and run the following command:
```bash
python 2_quantify.py
```

## Normalization

We usually do the CLR normalization for the cell by gene matrix before downstream analysis. This can be done in Seurat, here also provide a python script  `3_normalize.py` to do the CLR normalization.


