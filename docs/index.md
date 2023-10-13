# Covidscope reimplementation 

theme:
  name: readthedocs
  highlightjs: true
  hljs_languages:
    - yaml
    - rust

The reimplementation process is straightforward. You only need three separate files to start.

## Input file:

* Metadata
* Count matrix
* UMAP coordinates

Example input files are provided for your reference. Input files needed to be in .csv format.
We assume that users have analyzed their data in one of the most popular software Seurat in R. You can output your three input files from your Seurat object by:
* Metadata: your own metadata file in .csv format i.e. meta.csv
* Count matrix: 
```sh
# https://stackoverflow.com/questions/4558277/write-a-sparse-matrix-to-a-csv-in-r 
colnames(dense_matrix) <- c(“gene_name”, “barcode”,”expression”)
write.csv(dense_matrix, file=”matrix.csv”, row.names=FALSE) 
```

* UMAP:
If you have run the UMAP step, otherwise please refer to Seurat UMAP method and run it first. 
```sh
umap_coord <- dplyr::as_tibble(data.frame(seurat_object@reductions$umap@cell.embeddings), rownames = "id")
colnames(umap_coord) <- c(”id“, “UMAP1”,”UMAP2”)
write.csv(umap_coord, file=”umap.csv”, row.names=FALSE)
```

## Input file content:
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

## Dependencies:
* MongoDB installed
* Python installed
* MongoDB Compass installed

## Reimplementation steps:
1. Create a database named cov19atlas_new, and create three collections namely under the database:
* single_cell_meta_v4 
* umap
* matrix
2. Import the three datasets into MongoDB using MongoDB Compass
* single_cell_meta_v4 (meta.csv)
* umap (umap.csv)
* matrix (matrix.csv)

All data field during import needs to be in STRING format.

## Activate the virtual environment if you created

```sh
$ git clone ...
$ pip install -r requirements.txt
$ export FLASK_ENV=development
$ export FLASK_APP=manage.py
$ flask run
```
You will have local version of Covidscope running at 127.0.0.1:5000 by default.





