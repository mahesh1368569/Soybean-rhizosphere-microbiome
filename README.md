# Unraveling the microbiome dynamics in soybean rhizosphere region due to nitrogen and sulfur fertilization

This repository contains data and analysis workflows for a soybean experiment evaluating how nitrogen and sulfur treatments influence plant physiology, soil functional properties, and root/soil-associated microbial communities. The project integrates plant growth and physiology measurements with 16S rRNA bacterial/archaeal and ITS/ITS2 fungal microbiome datasets.
## Project overview

The main goal of this project is to connect soybean plant performance with soil biochemical indicators and microbiome responses under different nitrogen and sulfur treatment combinations.
- **64 soybean samples**
- **Nitrogen levels:** 0, 100, 150, and 200 (lbs/ac)
- **Sulfur levels:** 0, 15, 25, and 35 (lbs/ac)

The microbiome component includes:

- **16S rRNA amplicon data** for bacterial and archaeal community profiling
- **ITS/ITS2 amplicon data** for fungal community profiling
The workbook contains one sheet with 64 samples and 32 columns.

### Sample and treatment columns

| Column | Description |
|---|---|
| `Potnumber` | Pot or experimental unit number |
| `Sample ID` | Sample identifier used to link plant, soil, and microbiome data |
| `Exp_ID` | Experimental treatment or pot-level identifier |
| `Rep` | Replication/block identifier |
| `Nitrogen` | Nitrogen treatment level |
| `Sulphur` | Sulfur treatment level |

### Plant growth and yield-related variables

| Column | Description |
|---|---|
| `Plant Height` | Soybean plant height |
| `No. of Nodes` | Number of plant nodes |
| `No. of Pods` | Number of pods |
| `leaf area` | Leaf area measurement |
| `Pods dry weight` | Dry weight of pods |
| `Biomass (dry weight)` | Total aboveground dry biomass |
| `Trifolate Leaf weight (mg)` | Trifoliate leaf weight |

### Leaf physiology and spectral variables

| Column | Description |
|---|---|
| `E` | Transpiration rate |
| `A` | Photosynthetic assimilation rate |
| `gsw` | Stomatal conductance to water vapor |
| `ETR` | Electron transport rate |
| `Fs` | Steady-state fluorescence |
| `Fm'` | Maximum fluorescence under light-adapted conditions |
| `PhiPS2` | Effective quantum yield of photosystem II |
| `NDVI` | Normalized Difference Vegetation Index |
| `ChlM` | Chlorophyll index |
| `FlvM` | Flavonol index |
| `AnthM` | Anthocyanin index |
| `NFI` | Nitrogen balance or nitrogen-related fluorescence index |

### Soil functional and biochemical variables

| Column | Description |
|---|---|
| `AWCD` | Average well color development, commonly used as an indicator of microbial substrate utilization |
| `SDI` | Substrate diversity index |
| `POXC` | Permanganate oxidizable carbon |
| `BG` | β-glucosidase enzyme activity |
| `NAG` | N-acetyl-β-glucosaminidase enzyme activity |
| `ARS` | Arylsulfatase enzyme activity |
| `pH` | Soil pH |
