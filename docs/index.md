# Covidscope software architecture reimplementation guide 

The reimplementation process is straightforward. You only need three separate files to start.

## Input file

* Metadata
* Count matrix
* UMAP coordinates

Example input files are provided for your reference. Input files needed to be in .csv format.
We assume that users have analyzed their data in one of the most popular software Seurat in R. You can output your three input files from your Seurat object by:
* Metadata: your own metadata file in .csv format i.e. meta.csv
* Count matrix: 
```sh
# convert the sparse count matrix into dense matrix
# https://stackoverflow.com/questions/4558277/write-a-sparse-matrix-to-a-csv-in-r 
sparse2triples <- function(m) {
  SM = summary(m)
  D1 = m@Dimnames[[1]][SM[,1]]
  D2 = m@Dimnames[[2]][SM[,2]]
  data.frame(row=D1, col=D2, x=m@x)
}

# get matrix file
dense_matrix <- sparse2triples(srt_obj[["RNA"]]@counts)

colnames(dense_matrix) <- c("gene_name", "barcode","expression")
write.csv(dense_matrix, file="matrix.csv", row.names=FALSE) 
```
For large-scale single-cell data we also provide utilities to convert the 10X single-cell RNA-seq data to MongoDB import collection .csv format at [here](https://github.com/hiyin/covid19_cell_atlas_portal_supp/tree/main/utilities).

* UMAP:
If you have run the UMAP step, otherwise please refer to Seurat UMAP method and run it first. 
```sh
umap_coord <- dplyr::as_tibble(data.frame(seurat_object@reductions$umap@cell.embeddings), rownames = "id")
colnames(umap_coord) <- c(”id“, “UMAP1”,”UMAP2”)
write.csv(umap_coord, file=”umap.csv”, row.names=FALSE)
```

## Input file format
<table>
 <td>**File name**

   </td>
   <td>**Collection name**

   </td>
   <td>**Columns**

   </td>
  <tr>
   <td>metadata.csv

   </td>
   <td>single_cell_meta_v4

   </td>
   <td>["id", "meta_age_category", "meta_sample_id2","meta_patient_id", "meta_dataset", "level2", "meta_severity", "meta_days_from_onset_of_symptoms", "meta_outcome", "meta_gender", "Country"]

   </td>
  </tr>
  <tr>
   <td>matrix.csv

   </td>
   <td>matrix

   </td>
   <td>["barcode", "gene_name","expression"]

   </td>
  </tr>
  <tr>
   <td>umap.csv

   </td>
   <td>umap

   </td>
   <td>["UMAP1","UMAP2", "id"]

   </td>
  </tr>
</table>

## Database collection 

<table>
  <tr>
   <td>**Column name**

   </td>
   <td>**Value format**

   </td>
   <td>**Description**

   </td>
   <td>**Example**

   </td>
  </tr>
  <tr>
   <td colspan="4" >Collection schema: **single_cell_meta_v4**
   </td>
  </tr> 
   <td>id

   </td>
   <td>String

   </td>
   <td>Cell barcode (unique)

   </td>
   <td>
   </td>
  <tr>
   <td>meta_patient_id

   </td>
   <td>String

   </td>
   <td>Patient identifier

   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>meta_sample_id2

   </td>
   <td>String

   </td>
   <td>Sample identifier (one patient may have multiple samples)

   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>Country

   </td>
   <td>String

   </td>
   <td>Country origin of the dataset

   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>meta_age_category

   </td>
   <td>String

   </td>
   <td>Age category in intervals 

   </td>
   <td>e.g. “18-30”

   </td>
  </tr>
  <tr>
   <td>level2

   </td>
   <td>String

   </td>
   <td>Cell type predictions

   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>meta_severity

   </td>
   <td>String

   </td>
   <td>The severity of the symptoms regarding patient

   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>meta_dataset

   </td>
   <td>String

   </td>
   <td>Dataset origin

   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>meta_days_from_onset_of_symptoms

   </td>
   <td>Number

   </td>
   <td>Days from the onset of symptoms

   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>meta_gender

   </td>
   <td>String

   </td>
   <td>Gender 

   </td>
   <td>i.e. “female” or “male”

   </td>
  </tr>
  <tr>
   <td>meta_outcome

   </td>
   <td>String

   </td>
   <td>Health outcome

   </td>
   <td> i.e. “diseased”, “discharged” etc.

   </td>
  </tr>
   <td colspan="3" >Collection schema: **umap**

   </td>
   <td>
   </td>
  <tr>
   <td>id

   </td>
   <td>String

   </td>
   <td>Cell barcode (unique)

   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>UMAP1

   </td>
   <td>String

   </td>
   <td>X-coordinate of the cell

   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td>UMAP2

   </td>
   <td>String

   </td>
   <td>Y-coordiante of the cell

   </td>
   <td>
   </td>
  </tr>
  <tr>
   <td colspan="4" >Collection schema: **matrix**

   </td>
   </tr>
   <td>barcode

   </td>
   <td>String

   </td>
   <td>Cell barcode (unique)

   </td>
   <td>
   </td>
  <tr>
   <td>gene_name

   </td>
   <td>String

   </td>
   <td>Name of the gene expressed in the cell

   </td>
   <td>e.g. CD19

   </td>
  </tr>
  <tr>
   <td>expression

   </td>
   <td>Number

   </td>
   <td>
   </td>
   <td>
   </td>
  </tr>

</table>

## Dependencies
* MongoDB installed
* Python3.9 installed
* MongoDB Compass installed (optional, for visual operation of the database)
* sqlite installed

## Reimplementation one-by-one steps 
0. Install the system dependencies listed above, MongoDB, Python, MongoDB Compass (optional) and sqlite.
1. Create a database in MongoDB named cov19atlas_new, and create three collections namely under the database
    1. single_cell_meta_v4  
    2. umap                
    3. matrix              

2. Import the three datasets into MongoDB using MongoDB Compass, make sure all data field during import needs to be in STRING format.: 
    1. single_cell_meta_v4 ([importdata_meta.csv](https://covidscope-public-repository.s3.ap-east-1.amazonaws.com/user-tutorial/importdata_metadata.csv))
    2. umap ([importdata_umap.csv](https://covidscope-public-repository.s3.ap-east-1.amazonaws.com/user-tutorial/importdata_umap.csv)) and 
    3. matrix ([importdata_matrix.csv](https://covidscope-public-repository.s3.ap-east-1.amazonaws.com/user-tutorial/importdata_matrix.csv)).

3. Download Covidscope web portal resources below into flask_resources directory under your $HOME, without them the web portal can't be initialized
    1. [features.tsv](https://covidscope-public-repository.s3.ap-east-1.amazonaws.com/user-tutorial/features.tsv)
    2. [db.sqlite](https://covidscope-public-repository.s3.ap-east-1.amazonaws.com/user-tutorial/db.sqlite)

4. Clone the repository, install the packages for Covidscope and run the web portal code

## Quickstart guide (after data imported to MongoDB)
```sh
$ git clone https://github.com/hiyin/covid19_cell_atlas_portal.git
$ cd covid19_cell_atlas_portal

# create virtual environment
$ python3.9 -m pip install --upgrade pip
$ python3.9 -m venv venv # this installs the venv folder in the current directory
$ source venv/bin/activate
$ pip install -r requirements.txt

# deactivate the environment and reactivate it for a fresh start
$ deactivate
$ source venv/bin/activate

# initialize database download it and put it in flask_resources/ directory in yor $HOME
mv db.sqlite $HOME/flask_resources/
mv features.tsv $HOME/flask_resources/

# set environment
$ export FLASK_ENV=development
$ export FLASK_APP=manage.py

# launch
$ flask run
```
You will have local version of Covidscope running at 127.0.0.1:5000 by default.

You shall be able to have your own version of Covidscope running with your custom files if you could modify to the data format same as the example files we used here for demonstration!

## Prepare your data for Covidscope reimplementation
We assume that you start with the two common files after you have collected your single-cell RNA-seq data i.e. meta data, and a count matrix folder in 10X single-cell sequencing format. 

Below is a example pipeline to help you to process the files.
We asssume that you have a metadata file following our structures (if not please edit according to our example metadata)

1. Example metadata:
    [metadata](https://covidscope-public-repository.s3.ap-east-1.amazonaws.com/user-tutorial/importdata_metadata.csv)

2. Example 10X format matrix folder files:
    1. [barcodes.tsv.gz](https://covidscope-public-repository.s3.ap-east-1.amazonaws.com/raw/hoehn/barcodes.tsv.gz)
    2. [genes.tsv.gz](https://covidscope-public-repository.s3.ap-east-1.amazonaws.com/raw/hoehn/genes.tsv.gz)
    3. [matrix.mtx.gz](https://covidscope-public-repository.s3.ap-east-1.amazonaws.com/raw/hoehn/matrix.mtx.gz)
    
Download them and save as into a new directory named hoehn_2021/
```sh
meta <- read.csv("importdata_metadata.csv")
# Prepare matrix db collection source
hoehn <- Read10X("hoehn_2021/", gene.column = 1) # default tsv.gz files (downloaded from Covidscope)
srt_obj <- CreateSeuratObject(hoehn)
rownames(meta) <- meta$id
metadata_frame <- srt_obj@meta.data
univ_meta <- cbind(metadata_frame, meta)
srt_obj <- AddMetaData(srt_obj, univ_meta)

# if the user didn't have pca
srt_obj <- NormalizeData(srt_obj)
srt_obj <- FindVariableFeatures(srt_obj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(srt_obj)
srt_obj <- ScaleData(srt_obj, features=all.genes)
srt_obj <- RunPCA(srt_obj, features = VariableFeatures(object = srt_obj))
srt_obj <- RunUMAP(srt_obj, dims = 1:40)

# Prepare UMAP file for database import 
umap_coord <- dplyr::as_tibble(data.frame(srt_obj@reductions$umap@cell.embeddings
), rownames = "id")
colnames(umap_coord) <- c("id", "UMAP1","UMAP2")
write.csv(umap_coord, file="importdata_umap.csv", row.names = FALSE)

# Prepare matrix file for database import 
# input: a sparse matrix with named rows and columns (dimnames)
# returns: a data frame representing triplets (r, c, x) suitable for writing to a CSV file
sparse2triples <- function(m) {
  SM = summary(m)
  D1 = m@Dimnames[[1]][SM[,1]]
  D2 = m@Dimnames[[2]][SM[,2]]
  data.frame(row=D1, col=D2, x=m@x)
}

# get matrix file
towrite_matrix <- sparse2triples(srt_obj[["RNA"]]@counts)
colnames(towrite_matrix) <- c("gene_name", "barcode", "expression")
# write matrix data
write.csv(towrite_matrix, file="importdata_matrix.csv", row.names=FALSE) 
```

For a quick reproduction of Covidscope local version, you could download our prepared and processed files and directly import them into the MongoDB:
1. [matrix](https://covidscope-public-repository.s3.ap-east-1.amazonaws.com/user-tutorial/importdata_matrix.csv)
2. [metadata](https://covidscope-public-repository.s3.ap-east-1.amazonaws.com/user-tutorial/importdata_metadata.csv)
3. [umap](https://covidscope-public-repository.s3.ap-east-1.amazonaws.com/user-tutorial/importdata_umap.csv)

# How to use Covidscope data
For any user would like to use Covidscope data, downloaded from https://covidsc.d24h.hk/data, Covidscope provide the data in 10X format therefore you can read the three .gz files using Seurat package's function ```Read10X```. An example is given below
```sh
data <- Read10X("lee_2020/", gene.column = 1) # default tsv.gz files (downloaded from Covidscope)
srt_obj <- CreateSeuratObject(data)
```
