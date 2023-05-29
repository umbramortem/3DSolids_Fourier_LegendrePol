# 3DSolids_Fourier_LegendrePol
Desafortunadamente la presente plataforma (https://github.com/) no permite subi archivos muy grandes o pesados
mayores a 25mb, lo mismo que Optica en su plataforma y asi mismo los Journals de esta como Applied Opics. Por
ello los demas archivos mencionados en el manuscrito titulado "Tomography Fourier analysis for the construction 
of a 3D solid model for SARS-CoV-2 viral particles" como "HoWu_3Dmodel.mat" y "Ullah_3Dmodel.mat" se montaron 
en otra plataforma con la finalidad de estar al acceso de los lectores. Donde dicha plataforma solo depende de 
nuestra gestion, en google drive pero de gran capacidad. 

Ofrecemos una disculpa por tener que re-direccionar pero por cuestiones de economicas no cubrimos los costos de
GitHub para poder subir archivos de 1Gb de tamaño. Acontinuacion dejamos los siguientes enlaces de los archivos
mencionados:

1.- HoWu_3Dmodel.mat   link: https://drive.google.com/file/d/1M5LX7Jr3UiAIMbOHzR5IBOT-Q-x9RQpR/view?usp=sharing

2.- Ullah_3Dmodel.mat  link: https://drive.google.com/file/d/1pVUMo2qE1SbDKlK6fEISmIdODslu256P/view?usp=sharing

3.- Frames             link: https://drive.google.com/drive/folders/1OBWOUDH1ZaPz0YLEjntfqrOWLzhNx01o?usp=sharing


El ultimo enlace (carpeta "Frames") continee multiples imagenes y archivos de matlab (*.mat) para generar otros 
elementos y ejecutar sin problemas el codigo de "Results.m". Asi en "Frames" tenemos 4 carpetas, 2 para los 
archivos de Ho y 2 para Ullah. En las subcarpetas 40 Ho y 44 Ullah tenemos los frames (directos del video 1 y 2, 
respectivamente) en la carpeta frames, las mascaras binarias (carpeta Masks) y los segmentados finos (carpeta 
SRAS_Cell). Estas son necesarias en la linea 60 (Ho y 63 en el caso de los archivos de Ullah), 71 (74 para Ullah),
82 (85 Ullah) para el caso de los archivos en la subcarpeta frames. Luego los archivos de la carpeta Masks son 
necesarios para Ho en la linea 94 (97 para el caso de los archivos de Ullah).
Los archivos en la carpeta SARS_Cell son los necesarios para correr el codigo que genera los M imagenes del 
dolido 3D, nuestra propuesta en Python 2DLegendre_Fit. Asi, al contar con otra plataforma capaz de visualizar 
Solidos 3D importables, se tiene la posibilidad de obcervar nuestro resutado y no solo como lo proponemos en el 
manuscrito con la funcion "volumeViewer" de Matlab.


