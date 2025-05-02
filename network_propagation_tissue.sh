#!/bin/bash
#SBATCH --job-name tissue_network_prop
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --nodes 1
#SBATCH -a 1-38
#SBATCH -c 4
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu 8G
#SBATCH -o /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/tissue_network_prop-%j.o
#SBATCH -e /tscc/nfs/home/bsleger/bsl/SUD_cross_species/job_run_out/tissue_network_prop-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bsleger@ucsd.edu
#SBATCH --account csd795

source activate env-std-py38

cd /tscc/projects/ps-palmer/brittany/SUD_cross_species/tissue_networks

#brain_nervous_system=("amygdala" "basal_ganglion" "brain" "caudate_nucleus" "caudate_putamen" "central_nervous_system" "cerebellar_cortex" "cerebellum" "cerebral_cortex" "corpus_callosum" "corpus_striatum" "dentate_gyrus" "diencephalon" "forebrain" "frontal_lobe" "glia" "hippocampus" "hypophysis" "hypothalamus" "locus_ceruleus" "medulla_oblongata" "midbrain" "nervous_system" "neuron" "nucleus_accumbens" "occipital_lobe" "occipital_pole" "pons" "peripheral_nervous_system" "spinal_cord" "substantia_nigra" "subthalamic_nucleus" "telencephalon" "temporal_lobe" "thalamus" )
#adipose_tissue=("adipose_tissue")
#adrenal_gland=("adrenal_cortex" "adrenal_gland")
#blood_immune_system=("B_lymphocyte" "basophil" "blood" "blood_plasma" "blood_platelet" "dendritic_cell" "eosinophil" "granulocyte" "hematopoietic_stem_cell" "leukocyte" "lymphocyte" "macrophage" "mast_cell" "megakaryocyte" "monocyte" "mononuclear_phagocyte" "natural_killer_cell" "neutrophil" "serum" "T_lymphocyte" "thymocyte")
#bone_cartilage=("bone" "bone_marrow" "cartilage" "chondrocyte" "osteoblast" "tooth")
#cardiovascular_system=("aorta" "artery" "blood_vessel" "cardiac_muscle" "heart" "vascular_endothelial_cell" "vascular_endothelium")
#digestive_system=("cecum" "colon" "duodenum" "esophagus" "gastrointestinal_tract" "hepatocyte" "intestine" "ileum" "jejunum" "large_intestine" "liver" "pancreas" "pancreatic_islet" "rectum" "salivary_gland" "small_intestine" "stomach" "vermiform_appendix")
#endocrine_system=("corpus_luteum" "ovary" "pancreatic_islet" "placenta" "prostate_gland" "testis" "thyroid_gland")
#eye_ear=("choroid" "cochlea" "cornea" "ear" "eye" "lens" "retina")
#lymphatic_system=("lymph_node" "spleen" "thymus" "tonsil")
#muscular_system=("muscle" "skeletal_muscle" "smooth_muscle" "myometrium")
#reproductive_system=("mammary_epithelium" "mammary_gland" "ovarian_follicle" "ovary" "oviduct" "spermatid" "spermatocyte" "spermatogonium" "testis" "uterine_cervix" "uterine_endometrium" "uterus")
respiratory_system=("bronchial_epithelial_cell" "bronchus" "lung" "trachea" "global")
#skin_related=("epidermis" "hair_follicle" "keratinocyte" "skin" "skin_fibroblast")
#urinary_system=("kidney" "nephron" "podocyte" "renal_glomerulus" "renal_tubule" "urinary_bladder")
#other=("culture_condition_CD8+_cell" "embryo" "fetus" "global" "tear_gland" "trophoblast" "umbilical_cord" "umbilical_vein_endothelial_cell" "uroepithelium")


a=( "${brain_nervous_system[@]}" "${respiratory_system[@]}" )
echo ${#a[@]}

t=${a[$SLURM_ARRAY_TASK_ID-1]}


python  ../scripts/loco_network_prop.py $t
