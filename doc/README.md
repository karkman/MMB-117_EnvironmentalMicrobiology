# Wrangling of metadata

The original metadata was slightly modified. 

```{r}
metadata <- read_delim(
	"MMB117_Sample_sheet_2025.txt",
	delim = "\t",
	name_repair = 'minimal'
) %>%
	separate(SnowDepth, " ", into = "SnowDepth") %>%
	select(-ResponsiblePerson) %>% mutate(Sample = str_replace(Sample, "_", "-"))

write_delim(metadata, "doc/metadata.txt", delim="\t")
```