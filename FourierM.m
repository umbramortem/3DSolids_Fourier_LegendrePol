
function [Fou] = FourierM(Image, Mask)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This code constitutes an important tool for obtaining the results in
%    the manuscript entitled  "3D solid of SARS-CoV-2 viral particles 
%    applying Legendre polynomials from Tomography Fourier analysis", 
%    submitted to Journal of the Optical Society of America A (JOSA A) by
%    OOPTICA.
%
%    Correspondings Authors:
%    Dr. Jesus Alonso Arriaga Hernandez
%    jesus.arriagahdz@correo.buap.mx;    dr.j.a.arriaga.hernandez@gmail.com
%    Dr. Bolivia Teresa Cuevas Otahola
%    bolivia.cuevasotahola@viep.com.mx;          b.cuevas.otahola@gmail.com
%
%    This function contains three input data of square matrices of the form
%    MxM, which should be images transformed into Gray Tones (GT) of double
%    type.
%
%    In variable "X" we have the the medical image, in "Y" we have a binary
%    mask (SARS-CoV-2 cell). The output from the function (Fou) is the
%    unwrapped phase of the cell in the ROI of the medica image. We bear in 
%    mind that we aim to obtain the phase in all frames from the obtained 
%    phase of the information from a single frame in a 3D tomographic image,
%    to finally obtain a 3D model.
%
%    In this work, we analyze the majority of objects digitization or 3D
%    reconstruction techniques, considering from these the Fourier trnasfom
%    analysis by its simplicity. We induce a periodicity in the images,
%    filtering the periodic noise and finally unwrapping the phase. This
%    work was carried out in collaboration with Dr. Lilia Cedillo, Dr. Jorge
%    A Yañez Santos and Dr. Ygnacio Martínez, who greatly supported the
%    interpretation of the results. The videos of the medical image were
%    obtained in the research entitled: "Live imaging of SARS-CoV-2 
%    infection in mice reveals that neutralizing antibodies require Fc 
%    function for optimal efficacy";  with DOI: 10.1016/j.immuni.2021.08.015
%    and "Three-dimensional reconstruction by electron tomography for the
%    application to ultrastructural analysis of SARS-CoV-2 particles" with
%    DOI: https://doi.org/10.1007/s00795-021-00309-2  

%    ROI variables introduction (x) and binary mask of SARS-CoV-2 cell in
%    the ROI (y)

%    Depending on the input image type, if it is RGB unit8 we recommend
%    using the following lines, in case the function does not run properly.   

%imag = rgb2gray(Image);
%Image = double(imag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Optimal fringe size for an image represented by a square matrix (MxM)
%    Pattern that induces periodicity in the image to be studied.
%    We consider "N" as the fringes width in pixels, which in this case is
%    set to "N = 4" and its minimum value should be 1.

[A,B] = size(Image);
%    Fringe size in pixels in its optimal value (according to our case)
N = 4;

NN = 2*N;
n = A/NN;
aa = 2*floor(n);
patern = ones(A,B);

for i = 1 : aa

    a = 1 + ((i-1)*N);
    b = ((i)*N);

    if rem(i,2)==0

        patern(:, a:b ) = 0;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a, b] = size(Image);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z= Image.*patern;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Filter Construction and Location. Use of the "hanning" function to 
%    locate the first order of the two-dimensional Fourier Trasnform (FFT2)
%    in the proposed analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=size(z,1);
m=n;        
nh=n;
f=hanning(n);    ff=f;
[fx,fy]=meshgrid(f,ff);
s=fx.*fy;
m1=s;
fphi=z.*m1;
tf=fftshift(fft2(fphi));
[NFi,NCo]=size(tf);
C=zeros(NFi,NCo);
C(:,((NCo/2)+10):NCo)=tf(:,((NCo/2)+10):NCo);
B=C';
[vm,fil]=max(max(B));
[vm,col]=max(max(C));
[x,y]=meshgrid(-((n/2)-1):(n/2),-((m/2)-1):(m/2));
x0=col-(a/2);
y0=fil-(b/2);

%    Filter width
sigma=14;          

%    Filters, "filter_0" use for periodic noise and "filter" used in our
%    Fourier analysis
filter = exp(-1*((x-x0).^2+(y-y0).^2)/(2*sigma^2));
filter_0 = exp(-1*((x).^2+(y).^2)/(2*(1.5*sigma)^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Imag_test = Image;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    2D Fourier Tranform ordered towards the center
fft_TestIma = fftshift(fft2(Imag_test));

%    2D simple filter to remove periodic noise in the image ROI (Region of
%    Interest) in medical image (X) applying the 2D Fourier transform.
FT_fil = (fft_TestIma).*filter_0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XY_ifftimag = (ifft2(FT_fil));
XY_IfftImag = abs(XY_ifftimag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Removal and identification of the periodic noise
A = max(max(XY_IfftImag));
A2 = max(max(Imag_test));
if A < A2
    ruido = XY_IfftImag - double(Imag_test);
else
    ruido = double(Imag_test) - XY_IfftImag;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%                  The Fourier Transform Analysis
%%%%%%                         Our proposal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Periodicity induction with optimal fringes pattern

z2 = patern;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imagtest = (imresize(XY_IfftImag, [a b]));
z = z2.*imagtest;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    2D Fourier transform ordered towards the center, for the image ROI with
%    induced periodicity by the variable "patern"
FT_period_image = fftshift(fft2(z));
FT_patron = fftshift(fft2(z2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Image filtering by means of Fourier Analysis, for the ROI image with
%    periodicity

FT_fil_image = (FT_period_image).*filter;
FT_fil_patron = (FT_patron).*filter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XY_ift_Filimag = (ifft2(FT_fil_image));
XY_ift_patron = z2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Takeda Shifted
for i = 1 : a
    for j = 1 : b
        if z2(i,j) == 1
            Data_shiff(i,j) = 0;
        else
            Data_shiff(i,j) = 1;
        end
    end
end

FT_shif_patron = fftshift(fft2(Data_shiff));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FT_fil_image = (FT_period_image).*filter;
FT_fil_shif_patron = (FT_shif_patron).*filter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XY_ift_shif_patron = Data_shiff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Phase acquisition by means of the Fourier Analysis

S_A = ((XY_ift_shif_patron).*(XY_ift_Filimag));
S_B = ((XY_ift_patron).*(XY_ift_Filimag));
Super_product = (S_A)-(S_B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Result_0 = Mask.*abs(Super_product);

%    The following line is optional, since it could be removed, given that
%    it is intended only to show the resulting image.  Moreover, the graphic
%    type can be replaced with "mesh();", instead of the "imagesc();"
figure; imagesc(Result_0);
title('First final result with NAN');

%    Final result and output variable.
Fou = Result_0;

end
