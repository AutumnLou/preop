###########################################################################
# PhD
# 03_preop_therapy_graphics.R
###########################################################################

# Load dataset
source("scripts/02_data_comb_preop.R")

neoadj_hormone <- c("PCA0050", "PCA0103", "PCA0119", "PCA0176", "PCA0208")

# No Therapy - for order
order_data <- comb_data %>%
  filter(!(row.names(comb_data) %in% neoadj_hormone)) %>%
  select(ABCC11, ARMC4, CD38, DNAH8, DUSP18, EFCAB4B, ERO1LB, ESM1,
         FAM122C, FAP, FZD5, GRM8, HELB, NEU3, PI15, PLA2R1, SLC16A14, SOCS2,
         SPAG1, XBP1, ZHX3) %>%
  tibble::rownames_to_column(var = "patient")
  
median_ordered <- names(sort(robustbase::colMedians(as.matrix(order_data[,-1]))))

# All data
plot_data <- comb_data %>%
  #filter(!(row.names(comb_data) %in% neoadj_hormone)) %>%
  select(ABCC11, ARMC4, CD38, DNAH8, DUSP18, EFCAB4B, ERO1LB, ESM1,
         FAM122C, FAP, FZD5, GRM8, HELB, NEU3, PI15, PLA2R1, SLC16A14, SOCS2,
         SPAG1, XBP1, ZHX3) %>%
  tibble::rownames_to_column(var = "patient") %>%
  mutate(therapy = ifelse(patient %in% neoadj_hormone, "YES", "NO")) %>%
  arrange(therapy)


long_plot_data <- pivot_longer(plot_data, -c(patient, therapy), names_to = c("variable"))

long_therapy_data <- long_plot_data %>% filter(therapy == "YES")
  
ggplot(long_plot_data,
       aes(x=factor(variable, levels = median_ordered), y=value, group = patient)) +
  geom_line(aes(color = therapy)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Variable ()", y = "Value ()") +
  scale_color_manual(values=c('grey','black'))
  ggtitle("Figure 1")

## NOTE to self: Put the hormone patients at the end of the list!!
  ggplot() +
    geom_line(aes(x=factor(variable, levels = median_ordered), y=value, group = patient), 
              data = long_plot_data, colour = alpha("grey", 0.7)) +
    geom_line(size = 1.25, aes(x=factor(variable, levels = median_ordered), y=value, group = patient), 
              data = long_therapy_data) +  
    theme_bw() +
    labs(x = "mRNA Variables", y = "Normalized Transcript Expression") +
    theme(legend.position = "none", text = element_text(size=20),
          axis.text.x = element_text(angle=45, hjust=1))
  