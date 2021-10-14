# ress.py
Este código converte arquivos em coordenadas cartesianas para coordenadas cristalográficas. O arquivo de entrada deve ser o .out de um cálculo vc-relax das versões 6x do QuantumEspresso. 

As coordenadas do arquivo original devem estar em angstrom, e a extensão do novo arquivo deve ser ".res". Para torná-lo executável dê o comando, no terminal:

                                      chmod u+x ress.py
Criando um link simbólico em /usr/bin:

                                      sudo ln -sf $PWD/ress.py /usr/bin/ress.py

USO:
                                
                                      ress.py nome-do-arquivo.out
