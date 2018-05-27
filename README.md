# ESEM-ReplPack


===========================================
RUN the Java code as follows: 

java -Xmx5000m -jar WDP2.jar -dp "Your/DataPacks/Path/" -sp "Your/Save/Path/" -pc 13 | tee "Console-Output-to-File".txt"


============================================
Modified Weka Version 3.8.1

(AbstractInstance, Instance, SparseInstance, DenseInstance) in this version are modified in order to add additional fields to them. Specifically, an \textit{ID} field (of type Java String), an \textit{extra} field (of type Java String) and an InstIndex field (of type Java int) are added to the class AbstractInstance, the abstract Getter and Setter functions are added to class Instance, and these Getter and Getter functions are implemented in SparseInstance and DenseInstance classes
