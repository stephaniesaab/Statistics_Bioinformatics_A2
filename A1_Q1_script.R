# Load dataset
load("geneexpression2.rda")

# Inspect object names
ls()

# Rename object to "a1data"
a1data <- dat

# View structure
str(a1data)

# Dimensions
dim(a1data)

# Check first rows
head(a1data)

# Check row names
head(rownames(a1data))

# Check for missing values
sum(is.na(a1data))

# Summary statistics
summary(a1data)

# Extract row names
sample_names <- rownames(a1data)

# Split row names
split_names <- strsplit(sample_names, "_")

# Extract components
subject_status <- sapply(split_names, "[", 1)
cell_type <- sapply(split_names, "[", 2)

# Create metadata dataframe
metadata <- data.frame(
  Subject = subject_status,
  CellType = cell_type
)

head(metadata)
