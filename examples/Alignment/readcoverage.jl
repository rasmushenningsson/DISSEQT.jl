using SynapseClient
using DISSEQT.AlignUtils
using DISSEQT.SynapseTools

# Name of the run. (Should corrsespond to a subfolder of bamPath if in Synapse.)
runName = "H03UAAFXX"

# set to true to upload files
doUpload = true

# Should point to MyProject/Analysis/Alignment
alignmentFolder = "syn18694207" # VignuzziLabPublic/Projects/FitnessLandscapes/Analysis/Alignment


syn = SynapseClient.login()
coverageplots(syn, alignmentFolder, runName, doUpload, @__FILE__)


