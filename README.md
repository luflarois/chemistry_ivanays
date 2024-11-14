# chemistry_ivanays
Um pacote de teste para o solver - Migração para IA

Esse pacote lê um arquivo binário contendo dados de entrada do solver (rodas3_dyndt) da química do BRAMS.
O arquivo de entrada lido tem o seguinte nome composto:
**chem_inp_000030-0001.bin**
Após a leitura o programa faz o processo do solver e escreve um arquivo de saída com nome composto
**chem_std_000030-0001.bin**

|onde
|
|000030 = tempo da escrita (30.0 seg)
|0001 = Número do processador que foi usado

O arquivo é gerado por uma versão específica do BRAMS adaptada para esse fim. Esta versão se encontra em 
**https://github.com/luflarois/brams** no **branch brams_ai_chem**

O programa principal 


