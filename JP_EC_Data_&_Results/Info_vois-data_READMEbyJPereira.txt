10 VOIS ('lpmc','rpmc','lM1', 'rM1', 'lSMA', 'rSMA' <> 'visual1', 'broca1' _._ 'noise1', 'noise2');
	> ordenadas por import�ncia,
	> a seguir a '_._' s�o VOIs sem interesse fisiol�gico
	> a seguir a '<>' s�o VOIs com algum interesse fisiol�gico
10 Subjects;

2 Condi��es  - 'Decrease', 'Increase'

a) p/sujeito             -  8 blocos por run ; 12 pontos por bloco
b) p/ todos os sujeitos  -  80 blocos por run ('Decrease' + 'Increase') ; 12 pontos por bloco

Caso b):
1 - tarefa 1 NF_BEF_DEC   -  VOIS antes de NF  -  1run  -  960 time points (  80 blocos 'Decrease')
2 - tarefa 1 NF_BEF_INC   -  VOIS antes de NF  -  1run  -  960 time points (  80 blocos 'Increase') 

5 - tarefa 2 NF_DEC       -  VOIS durante NF   -  3runs - 2880 time points ( 240 blocos 'Decrease')
6 - tarefa 2 NF_INC       -  VOIS durante NF   -  3runs - 2880 time points ( 240 blocos 'Increase')

3 - tarefa 3 NF_AFTER_DEC -  VOIS depois de NF -  1run  -  960 time points (  80 blocos 'Decrease')
4 - tarefa 3 NF_AFTER_INC -  VOIS depois de NF -  1run  -  960 time points (  80 blocos 'Increase') 

------------------------------------------------------------
-> As quest�es seriam : (poss�veis hip�teses: numero de blocos x n�mero de pontos __ n�mero de VOIs)

TEST_1
comparar tarefa 2( Decrese) com tarefa 2 (Increase) - teste em n�o resting state, variabilidade inter-sujeito,
						      como controlar para este fen�meno?

	>>> hipotese A  -   1 x 2880     vs  1 x 2880     __ [3:6]
	>>> hipotese B  -  80 x   12     vs 80 x   12     __ [3:6]  (n�o fazer baseado no nosso artigo)
	>>> hipotese C  -   4 x  770     vs  4 x  770     __ [3:6]
	>>> hipotese D  -  25 x  108     vs 25 x  108     __ [3:6]
	>>> Teste Princ -  24 x   48(36) vs 24 x   48(36) __ [3:6]

TEST_2
comparar tarefa 1( Decrease) com tarefa 3( Decrease) - teste em resting state

	>>> hip�tese A  -  1 x 960  vs  1 x 960 __ [3:6]
	>>> hip�tese B  -  2 x 480  vs  2 x 480	__ [3:6]	
	>>> Teste Princ - 30 x  32  vs 30 x  32 __ [3:6]

TEST_3
comparar tarefa 1( Increase) com tarefa 3(Increase)

	>>> hip�tese A  -  1 x 960  vs  1 x 960 __ [3:6]
	>>> hip�tese B  -  2 x 480  vs  2 x 480	__ [3:6]	
	>>> Teste Princ - 30 x  32  vs 30 x  32 __ [3:6]	

TEST_4 possivelmente:
comparar Tarefa 2 ('Decrease' / 'Increase') com Tarefa 1+3 ('Decrease' / 'Increase')  [3runs vs 2runs]

	>>> hip�tese A  -  1 x 2880  vs 1 x 1920 __ [3:6]
	>>> hip�tese B  -  2 x 1440  vs 2 x  960 __ [3:6]
	>>> hip�tese C  -  4 x  770  vs 4 x  480 __ [3:6]
