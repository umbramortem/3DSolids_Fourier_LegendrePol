# 3DSolids_Fourier_LegendrePol

English Instructions

Unfortunately, this platform (https://github.com/) does not allow uploading very large or heavy files larger 
than 25mb, the same as Optica on its platform and likewise its Journals as Journal of the Optical Society of
America A (JOSA A) by OPTICA. Therefore, the other files mentioned in the manuscript entitled "3D solid of 
SARS-CoV-2 viral particles applying Legendre polynomials from Tomography Fourier analysis" such as 
"HoWu_3Dmodel.mat" and "Ullah_3Dmodel.mat" were mounted on another platform with the purpose of being 
accessible to readers. Where said platform only depends on our management, in google drive but with great 
capacity.

We apologize for having to redirect but for economic reasons we do not cover the costs of GitHub to be
able to upload files of 1Gb in size. Next we leave the following links of the files mentioned:

1.- HoWu_3Dmodel.mat   link: https://drive.google.com/file/d/1M5LX7Jr3UiAIMbOHzR5IBOT-Q-x9RQpR/view?usp=sharing

2.- Ullah_3Dmodel.mat  link: https://drive.google.com/file/d/1pVUMo2qE1SbDKlK6fEISmIdODslu256P/view?usp=sharing

3.- Frames             link: https://drive.google.com/drive/folders/1OBWOUDH1ZaPz0YLEjntfqrOWLzhNx01o?usp=sharing

The last link ("Frames" folder) contains multiple images and matlab files (*.mat) to generate others.
elements and execute without problems the code of "Results.m". So in "Frames" we have 4 folders, 2 for the
files for Ho and 2 for Ullah. In the subfolders 40 Ho and 44 Ullah we have the frames (direct from video 1 and 2,
respectively) in the frames folder, binary masks (Masks folder) and fine segments (
SARS_Cell). These are needed in line 60 (Ho and 63 in the case of Ullah files), 71 (74 for Ullah),
82 (85 Ullah) for the files in the frames subfolder. Then the files in the Masks folder are
necessary for Ho in line 94 (97 in the case of Ullah files).
The files in the SARS_Cell folder are the ones needed to run the code that generates the M images of the
pained 3D, our proposal in Python 2DLegendre_Fit. Thus, having another platform capable of visualizing
Importable 3D solids, it is possible to see our result and not only as we propose it in the
manuscript with the function "volumeViewer" of Matlab.
Finally, in the main folder "Frames" there are the last two folders "562 Hon Wu" and "660 Ullah"
contain the 562 fit files for Ho (660 in the case of Ullah) obtained with 2DLegendre_Fit. In this case
it only remains the case of graphing it. For this reason, the importance of attaching this data, since otherwise
it does not apply. The program because it is very specific to the cases, that is to say, that only the code 
"Results.m" is exclusive for this investigation in the case of the tomographic studies in Visualization 1 and 
Visualization 2. While the "FourierM.m" and "2DLegendre_Fit" codes are general and applicable to multiple cases, 
we only have to consider the illumination in the tomographic image, since it should not have sweeps or shadows 
due to the inclination of the illumination.


Instrucciones en español

Desafortunadamente la presente plataforma (https://github.com/) no permite subi archivos muy grandes o pesados
mayores a 25mb, lo mismo que Optica en su plataforma y asi mismo los Journals de esta como Optical Society of
America A (JOSA A) perteneciente a OPTICA. Por ello los demas archivos mencionados en el manuscrito titulado 
"3D solid of SARS-CoV-2 viral particles applying Legendre polynomials from Tomography Fourier analysis" como 
"HoWu_3Dmodel.mat" y "Ullah_3Dmodel.mat" se montaron en otra plataforma con la finalidad de tener el acceso de 
los lectores. Donde dicha plataforma solo depende de nuestra gestion, en google drive pero de gran capacidad. 

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
Finalmente, en la carpeta principal "Frames" se encuantran las ultimas dos carpetas "562 Hon Wu" y "660 Ullah" 
contienen los 562 archivos de ajuste para Ho (660 en el caso de Ullah) obtenidos con 2DLegendre_Fit. En este caso
solo resta el caso de graficarlo. Por ello la importancia de anexar estos datos ya que de lo contrario no corre 
el programa porque es muy especifico de los casos, es decir, que solo el codigo "Results.m" es exclusivo para 
esta invetsigacion en el caso de los estudios tomograficos en Visualization 1 y Visualization 2. Mientras que los 
codigos "FourierM.m" y "2DLegendre_Fit" son generales y aplicables s multiples casos, solo debemos considerar la
iluminacion en la imagen tomografica, ya que no debe tener barridos ni sombras por inclinacion de la iluminacion. 


