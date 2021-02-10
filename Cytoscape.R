#Install all relevant packages
######################################

if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}

if(!"igraph" %in% installed.packages()){
  install.packages("igraph")
}

if(!"plyr" %in% installed.packages()){
  install.packages("plyr")
}

if(!"ndexr" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("ndexr")
}

######################################

#Load in libraries
library(RCy3)
library(ndexr)
library(plyr)
library(igraph)

#You can start with NDEx by first establishing a connection.
ndexcon <- ndex_connect()

############################################################################
                          #  DEGs visualisation #
############################################################################

######################################
#OsRKD DEGs 
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD geneNet 0.15")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 


######################################
#TaRKD DEGs 
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD TaRKD DEGs")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 


######################################
#Microspore DEGs
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Microspore DEGs")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 


######################################
#Tillering DEGs
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Tillering DEGs")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 



######################################
#NGR5 DEGs
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD NGR5 DEGs")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 


############################################################################
                        #  geneNet visualisation #
############################################################################

######################################
#TaRKD geneNet 0.15 
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD TaRKD geneNet 0.15")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 


######################################
#Microspore geneNet 0.15 
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Microspore geneNet 0.15")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 


######################################
#Tillering geneNet 0.15 
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Tillering geneNet 0.15")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 



######################################
#NGR5 geneNet 0.15 
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD NGR5 geneNet 0.15")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 



######################################

############################################################################
                    #  Tissue expression visualisation #
############################################################################

######################################
#Pre-emergence inflorescence
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Pre-emergence inflorescence")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 




######################################
#Post-emergence inflorescence
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Post-emergence inflorescence")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 




######################################
#Anther
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Anther")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 




######################################
#Pistil
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Pistil")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 



######################################
#Embryo 25 DAP
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Embryo 25 DAP")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 




######################################
#Endosperm 25 DAP
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Endosperm 25 DAP")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 




######################################
#Seed 5 DAP
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Seed 5 DAP")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 




######################################
#Seed 10 DAP
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Seed 10 DAP")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 



######################################
#Leaves 20 days
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Leaves 20 days")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 





######################################
#Seedling 4 leaf stage
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Seedling 4 leaf stage")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 



######################################
#Shoots 25 DAP
######################################

#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Shoots 25 DAP")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 


#Search for public network OsRKD
networks <- ndex_find_networks(ndexcon, "OsRKD Seedling 4 leaf stage")
print(networks[,c("name","externalId","nodeCount","edgeCount")])

#Since only one has popped up we can go ahead to the get network function
networkId = networks$externalId[1]
network = ndex_get_network(ndexcon, networkId)
print(network)

#Importing network into R
cytoscapePing()
importNetworkFromNDEx(networkId) 
