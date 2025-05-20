# Soil Legacies of Genotypic Diversity Promote Plant Performance by Mediating Ontogenetic Allometry

This repository contains the code and data for the manuscript:  
**“Soil legacies of genotypic diversity promote plant performance by mediating ontogenetic allometry.”**

---

## 📁 Project Structure

- **`GeneticDiversity.Rproj`**  
  RStudio Project file. Open this to launch the working environment for the project.

---

## 📊 Data Files

- **`Soil.xlsx`**  
  Contains soil samples collected at the end of the **training phase**.  
  For each pot, soil was sampled from 12 randomly selected points using sterile medical syringes and then combined into a single composite sample.  
  This dataset is used to assess soil legacies resulting from different genotypic diversity treatments.

- **`DataSecond.xlsx`**  
  Contains plant trait data collected at the end of the **testing phase**.  
  Plants were cut at the root base; aboveground parts were harvested, and belowground parts were washed free of soil.  
  We recorded plant height and counted clonal ramets, inflorescences, seeds, and corms in each pot.  
  All plant materials were then oven-dried at 70°C for 72 hours to a constant weight.  
  We measured seed biomass, inflorescence biomass (including seed biomass), aboveground vegetative biomass, fine root biomass, corm biomass, and rhizome biomass.

- **`MonthlyData.xlsx`**  
  Contains monthly monitoring data of *S. mariqueter* phytometers from April 22 to October 22, 2022.  
  This dataset was used to assess plant growth and reproductive performance in soils trained by different genotypic diversity levels.

- **`PlotTreatment.xlsx`**  
  Metadata linking each pot ID to its corresponding experimental treatment.  
  Includes information such as genotypic diversity levels and treatment groupings.

---

## 🧪 Code Files

- **`GeneticDiversity.R`**  
  Main R script used in the manuscript.  
  Performs data import, cleaning, statistical analysis, and visualization.  
  This script generates the main results and figures presented in the paper.

- **`stepba function`**  
  A custom R function used for **model simplification** via stepwise backward selection.  
  Helpful in identifying minimal adequate models for statistical analysis.

---

## ▶️ How to Use

1. Open `GeneticDiversity.Rproj` in RStudio.
2. Ensure the data file is in the same directory as the R scripts.
3. Run `GeneticDiversity.R` to reproduce the analyses and generate figures used in the manuscript.

