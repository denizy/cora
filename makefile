CC=g++
CFLAGS=-static -march=native -O3 -o 

all: collapse collapse_L cora fastqSplit linkConstruct homTable_setup homTable_setup_L mappingInference mappingInference_L mappingInference_E3 mappingInference_L_E3 mappingInference_E4 mappingInference_L_E4 mappingInference_E5 mappingInference_L_E5 mappingInference_E6 mappingInference_L_E6

collapse: collapse.cpp
	$(CC) $(CFLAGS) collapse collapse.cpp
collapse_L: collapse.cpp
	$(CC) -DCHR_SHORT  $(CFLAGS) collapse_L collapse.cpp
cora: cora.cpp
	$(CC) $(CFLAGS) cora cora.cpp
fastqSplit: fastqSplit.cpp
	$(CC) $(CFLAGS) fastqSplit fastqSplit.cpp
linkConstruct: linkConstruct.cpp
	$(CC) $(CFLAGS) linkConstruct linkConstruct.cpp
homTable_setup: homTable_setup.cpp
	$(CC) $(CFLAGS) homTable_setup homTable_setup.cpp -fopenmp
homTable_setup_L: homTable_setup.cpp
	$(CC) -DCHR_SHORT $(CFLAGS) homTable_setup_L homTable_setup.cpp -fopenmp 
mappingInference: mappingInference.cpp
	$(CC) $(CFLAGS) mappingInference mappingInference.cpp
mappingInference_L: mappingInference.cpp
	$(CC) -DCHR_SHORT $(CFLAGS) mappingInference_L mappingInference.cpp
mappingInference_E3: mappingInference.cpp
	$(CC) -DED_THREE $(CFLAGS) mappingInference_E3 mappingInference.cpp
mappingInference_L_E3: mappingInference.cpp
	$(CC) -DCHR_SHORT -DED_THREE $(CFLAGS) mappingInference_L_E3 mappingInference.cpp
mappingInference_E4: mappingInference.cpp
	$(CC) -DED_FOUR $(CFLAGS) mappingInference_E4 mappingInference.cpp
mappingInference_L_E4: mappingInference.cpp
	$(CC) -DCHR_SHORT -DED_FOUR $(CFLAGS) mappingInference_L_E4 mappingInference.cpp
mappingInference_E5: mappingInference.cpp
	$(CC) -DED_FIVE $(CFLAGS) mappingInference_E5 mappingInference.cpp
mappingInference_L_E5: mappingInference.cpp
	$(CC) -DCHR_SHORT -DED_FIVE $(CFLAGS) mappingInference_L_E5 mappingInference.cpp
mappingInference_E6: mappingInference.cpp
	$(CC) -DED_SIX $(CFLAGS) mappingInference_E6 mappingInference.cpp
mappingInference_L_E6: mappingInference.cpp
	$(CC) -DCHR_SHORT -DED_SIX $(CFLAGS) mappingInference_L_E6 mappingInference.cpp
clean:
	rm collapse collapse_L cora fastqSplit linkConstruct homTable_setup homTable_setup_L mappingInference mappingInference_L mappingInference_E3 mappingInference_L_E3 mappingInference_E4 mappingInference_L_E4 mappingInference_E5 mappingInference_L_E5 mappingInference_E6 mappingInference_L_E6
