raw_data_labels = read.table("Label_abundance_IBD.csv", h = F, sep = ",")
raw_data_labels_vert = as.data.frame(t(raw_data_labels))
colnames(raw_data_labels_vert) = c("ech","labels_pre_SPARTA")
raw_data_labels_vert

split_labels = read.table("Annotation_samples_separation_Iteration_0_run2.csv", h = F, sep = ",")
split_labels_vert = as.data.frame(t(split_labels))
colnames(split_labels_vert) = c("ech","labels_post_SPARTA")
split_labels_vert

df = merge(raw_data_labels_vert, split_labels_vert, by = "ech")
df
