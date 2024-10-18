# Palettes

#### Lakes ####
palette_lake <- c("#00A651", "#8DC63F", "#47B648", "#ED1C24", "#00AEEF", "#2E3192")
names(palette_lake) <- c("L226N", "L226S", "L226", "L227", "L224", "L373")

palette_stratum <- c("#55D6BE", "#F08700", "#2E5EAA", "#2E5EAA")
names(palette_stratum) <- c("Epilimnion", "Metalimnion", "Hypolimnion", "Profundal")

palette_stratum_shape <- c(17, 15)
names(palette_stratum_shape) <- c("Epilimnion", "Metalimnion")

palette_manipulation <- c("#E2F0D9", "#A9D18E", "#C5E0B4", "#F4F9F0",
                          "#A9D18E", "#A9D18E", "#b9c4cb", "#5CD3FF")
names(palette_manipulation) <- c("P only", "13N:P", "5N:P", "Fish",
                                 "C+N+P", "C+N", "Metals", "Drawdown")


#### Time periods ####
# Seasons
palette_season <- c("#63D2FF", "#78BC61", "#FFE66D", "#FF6B6B")
names(palette_season) <- c("Winter", "Spring", "Summer", "Fall")

palette_month <- c("#63D2FF", "#99E2FF", "#C2EEFF",
                   "#78BC61", "#91C87E", "#A9D49B",
                   "#FFE66D", "#FFEE99", "#FFF5C2",
                   "#FF6B6B", "#FF9999", "#FFC2C2")
names(palette_month) <- c("January", "February", "March",
                          "April", "May", "June",
                          "July", "August", "September",
                          "October", "November", "December")

palette_stratification <- c("#1E90FF", "#FF8D1E")
names(palette_stratification) <- c("Stratified", "Mixed")

#### Taxonomy and functions ####
# Opisthokonts, land plants, and "protists"
palette_taxgroup <- c("#FF2D56", "#D81159", "#00B26E", "yellow")
names(palette_taxgroup) <- c("Metazoa", "Fungi", "Plants", "Other (protists)")

# Eukaryotic supergroups
palette_supergroup <- c("#75DAFF", "#FF2D56", "#3DFFE5", "#FF7CEB", "#00E288",
                        "#F6FF00", "#8A84E2", "#E8E8E8", "#CC7D5D")
names(palette_supergroup) <- c("Stramenopiles", "Opisthokonta", "Hacrobia", "Alveolata", "Archaeplastida",
                               "Rhizaria", "Excavata", "Amoebozoa", "Apusozoa")

# Eukaryotic subdivisions
palette_subdivision <- c("#75DAFF",
                         "#FF007B",
                         "#D81159",
                         "#FF7CEB",
                         "#F6FF00",
                         "#EA00C7",
                         "#00E288",
                         "#CE00B6",
                         "#276ED8",
                         "#FFDAF5",
                         
                         "#3DFFE5",
                         "#0099AA",
                         "#FF8787",
                         "#9BFFEB",
                         "#00B26E",
                         "#31E2DF",
                         "#8A84E2",
                         "#CC7D5D",
                         
                         "#767676",
                         "#A9A9A9",
                         
                         "black",
                         
                         "#EEEEEE")
names(palette_subdivision) <- c("Gyrista",
                                "Apicomplexa",
                                "Fungi",
                                "Dinoflagellata",
                                "Cercozoa",
                                "Ciliophora",
                                "Chlorophyta",
                                "Perkinsea",
                                "Bigyra",
                                "Chrompodellids",
                                
                                "Cryptophyta",
                                "Kathablepharida",
                                "Choanoflagellata",
                                "Centroplasthelida",
                                "Streptophyta",
                                "Haptophyta",
                                "Euglenozoa",
                                "Collodictyonidae",
                                
                                "Unclassified Opisthokonta",
                                "Unclassified Alveolata",
                                
                                "Unclassified subdivision",
                                
                                "Other subdivision")


# Phototroph classes
palette_phototroph <- c("#75DAFF",
                        "#0096CC",
                        "#00B4F5",
                        "#1FC3FF",
                        "#47CEFF",
                        "#99E4FF",
                        "#C2EFFF",
                        "#EBFAFF",
                        "#FF7CEB",
                        "#00E288",
                        "#C2FFE7",
                        "#8BE271",
                        "#D2FF91",
                        "#00B26E",
                        "#8A84E2")
names(palette_phototroph) <- c("Chrysophyceae",
                               "Coscinodiscophyceae",
                               "Bacillariophyceae",
                               "Phaeothamniophyceae",
                               "Bolidophyceae",
                               "Eustigmatophyceae",
                               "Dictyochophyceae",
                               "Xanthophyceae",
                               "Dinophyceae",
                               "Chlorophyceae",
                               "Trebouxiophyceae",
                               "Mamiellophyceae",
                               "Chlorodendrophyceae",
                               "Zygnemophyceae",
                               "Euglenida")

palette_phytoplankton <- c("grey85", "#00E288", "#75DAFF", "#FF7CEB", "#3DFFE5",
                           "#C9F7FF", "#8A84E2", "black")
names(palette_phytoplankton) <- c("Cyanobacteria", "Chlorophytes", "Chrysophytes", "Dinoflagellates", "Cryptophytes",
                                  "Diatoms", "Euglenophytes", "total_phytoplankton")


# Trophic functional groups
palette_function <- c("#EDD55C", "#AB6B42",
                      "#F55D3E", "#B78BD4",  # orange soda, lavender floral, 
                      "#76BED0", "#56E4C0", "#62BE81",  # dark sky blue, medium aquamarine, emerald
                      "#FDFF80"   # arylide yellow, laser lemon, brown sugar
)
names(palette_function) <- c("symbiotrophs", "pathotrophs",
                             "parasites", "heterotrophs/parasites",
                             "heterotrophs", "heterotrophs/mixotrophs", "mixotrophs",
                             "phototrophs")
