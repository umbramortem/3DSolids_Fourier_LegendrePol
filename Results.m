
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    This code constitutes an important tool for obtaining the results in
%    the manuscript entitled  "Tomography Fourier analysis for the
%    construction of a 3D model for SARS-CoV-2 viral particles Phase
%    analysis simulating the Takeda method to obtain a 3D profile of
%    SARS-CoV-2 viral particles", submitted to Applied Optics by OPTICA
%    (former OSA-Optical Society of America).
%
%    Recalling that the ROI are built from 3D tomographic image data
%    from the research work "Live imaging of SARS-CoV-2 infection in mice
%    reveals that neutralizing antibodies require Fc function for optimal
%    efficacy" with DOI: 10.1016/j.immuni.2021.08.015 and from the work 
%    entitled "Three-dimensional reconstruction by electron tomography for
%    the application to ultrastructural analysis of SARS-CoV-2 particles"
%    with DOI: 10.1007/s00795-021-00309-2
%    
%    From such work, we selected several frames (video image captures from
%    the tomographic image) and for each frame we obtained the phase of the
%    observed SARS-CoV-2 cells, segmented by simulating the Takeda Method
%    (fringes projection along with the use of the Fourier Transform).
%    Subsequently, we integrate the phase of the different SARS-CoV-2 cells
%    frames in a 3D matrix, where we obtain the 3D object of the SARS-CoV-2
%    virion particle.
%
%    Correspondings Authors:
%    Dr. Jesus Alonso Arriaga Hernandez
%    jesus.arriagahdz@correo.buap.mx;    dr.j.a.arriaga.hernandez@gmail.com
%    Dra. Bolivia Teresa Cuevas Otahola
%    bolivia.cuevasotahola@viep.com.mx;          b.cuevas.otahola@gmail.com
%
%    This algorithm contains a simple routine to call the main files 
%    required to obtain the results as well as to apply our main function
%    "FourierM.m".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    We should consider in the first place a period pattern (fringes) that
%    induce periodicity in the images. Such pattern is built in the function
%    "FourierM.m", which also takes into account the optimal fringes pattern
%    according to size of the images under study for the cells in the
%    studied videos.  The instruction "load()" works according to the file
%    location in the work PC.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Previously, all frames from the ROI and masks segmenting cells 1 and 2
%    are saved, in each work corresponding to the SARS-CoV-2 tomographic 
%    image obtained from the video of the tomographic images. Subsequently,
%    we open the image files of the ROI using the following routines 
%    (which opes all the images sorted with a single instruction at the 
%    beginning of the following loop).
%
%    Instructions to read and open the ROI of cells 1 and 2 and binary mask
%    corresponding to each frame, where the image and mask sizes should be
%    equal.
%
%    Cell 1 
Files = dir('C:\Location on your Computer\NEW WORK\Hong Wu\Frames\*.png');

%    Cell 2
%Files = dir('C:\Location on your Computer\NEW WORK\Ullah\Frames\*.png');

J = length(Files);

for k = 1:J
    %    Instruction to sort the frames, from the highest to the lowest 
    %    part of the SARS-CoV-2 cell.
    %    Cell 1 Frames 
    File_png = strcat('C:\Location on your Computer\NEW WORK\Hong Wu\Frames\Frame_', num2str(k), '.png');

    %    Cell 2 Frames
    %File_png = strcat('C:\Location on your Computer\NEW WORK\Ullah\Frames\Frame_', num2str(k), '.png');

    Image= imread(File_png);

    %    Instruction to segment with a square geometry the ROI, resulting 
    %    in a frame containing information from the desired cell only.

    %    Cell 1
    cropp_image = imcrop(Image, [1047.5 304.5 280 280]);

    %    Cell 2
    %cropp_image = imcrop(Image, [1112 411 219 219]);

    frames(:,:,:,k) = cropp_image;

    %    Instruction to read the binary masks segmenting according to the
    %    information up to the SARS-CoV-2 vision particles
    %    borders.

    %    Binary mask corresponding to cell 1
    mask = strcat('C:\Location on your Computer\NEW WORK\Hong Wu\Masks\Mask_', num2str(k), '.mat');

    %    Binary mask corresponding to cell 2
    %mask = strcat('C:\Location on your Computer\NEW WORK\Ullah\Masks\Mask_', num2str(k), '.mat');

    masks = load(mask);
    Masks(:,:,k) = cell2mat(struct2cell(masks));

end
[a,b,c] = size(frames(:,:,:,J));

%    Instruction for the medical image or frames in Cell 1, Hong Wu work
A = 2*a;                         B = 2*b;

%    Instruction for the medical image or frames in Cell 2, Ullah work
%A = 3*a;                         B = 3*b;

for k = 1:J
%   Instruction to increase the size of the ROI images to optimize the
%   3D model visualization
    frame = imresize(frames(:,:,:,k), [A B]);
        
    %    Results from the phase unwrapping of the SARS-CoV-2  Cell 1 
    %    information for each frame. These instructions are different since
    %    the function admits as input variable in X the images with their
    %    corresponding masks in double and GT format, whereas the variable Y
    %    can be constant. However, variable Z should always be the proposed
    %    pattern.e ser el patron propuesto
    
    imag = rgb2gray(frame);
    Image = double(imag);
    Aux = Image.*Masks(:,:,k);
    Result(:,:,k) = FourierM(Aux,1);
    ResultB(:,:,k) = Result(:,:,k).*Masks(:,:,k);

    %    In the following routine the edges with values equal to zero in the
    %    binary mask are removed, segmenting more precisely the cell
    %    information.  This constitutes the first version of the 3D model
    %    from the faithful data of the cell.
   
    Aux = ResultB(:,:,k);
    for i=1:A
        for j=1:B
            if Aux(i,j) == 0
                Result_0(i,j,k) = NaN;
            else
                Result_0(i,j,k) = Aux(i,j);
            end
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Routine to generate a preliminar version of the matrix with shape mxmxz
%    with the frames obtained from the Video. We build the elements from the
%    corresponding level curves or isophotes, integrating the data in 3D
%    matrix.
%    3D  Matrix 562x562x44, for cell 1 (Hong Wu work) and 660x660x40 for
%    cell 2 (Ullah work)

%    Rule of three for the scale, 50nm for cell 1 and 100nm for cell 2, in 
%    the tomographic image to resolve the measurement of the diameters
%    Cell 1
XX = (100/993)*943;

%    Cell 2
%XX = (100/1576)*1168;

[X, Y]=meshgrid(0:(XX/A):XX);
X=imresize(X,[A B]);    Y=imresize(Y,[A B]);

figure;
hold on
for ii = 1 : J
    zz = (3.5*J) - (3.5*ii);         %  separation between multiple Z planes 
    [~,h] = contour(X,Y,Result_0(:,:,ii),30);   % plot contour at the bottom
    h.ContourZLevel = zz;     
end
hold off
view(3);
xlabel('X [nm]'); ylabel('Y [nm]'); zlabel('Z [nm]', 'Rotation',0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    The latter routines are instructions to segment using the toolbox IMAGE
%    PROCESSING AND COMPUTER VISION, and segmenting with Image Segmenter in
%    MATLAB, in addition to the processing functions we are presenting, 
%    where we simulate the Takeda method to induce periodicity in medical
%    images and to recover their phase.
%    Subsequently, in the following section, we apply our propose routine
%    "2DLegendre_Fit" developed in Python 3, which uses previous results,
%    exported in .mat files in matlab to Python to carry out data 
%    interpolation using the 1D Legendre polynomials which generates 
%    neighboring images from the data recovered phase of SARS-CoV-2 cells.
%    Two phase images are considered to generate an intermediate one, which
%    is subsequently repeated until as many images as required are generated
%    or as the algorithm requires, as we mention in the manuscript. Where
%    "2DLegendre_Fit"is a modified algorithm from the first version
%    published in the paper "Two-dimensional Legendre polynomials as a basis
%    for interpolation of data to optimize the solution of the irradiance
%    transport equation analyzed as a boundary problem on surfaces testing"
%    with DOI: 10.1364/AO.58.005057. Since our main goal is to build a 3D
%    model of the SARS-CoV-2 virion particles, we build MxMxM matrices (in 
%    our case according the previous AxAxA data), specifically 562x562x562 
%    for cell 1 and 660x660x660 for cell 2. 
%    Generating 562 from 44 and 660 from 40 for cells 1 and 2, respectively.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    Instruction to read the interpolation images to generate the 3D object
%    from the 3D matrix: 
%    562 images for cell 1

read_FitImages = dir('C:\Location on your Computer\NEW WORK\Interpoladas\Hong Wu\*.mat'); 

%    660 images for cell 2
%read_FitImages = dir('C:\Location on your Computer\NEW WORK\Interpoladas\Ullah\*.mat');

M = length(read_FitImages);

for m = 1:M
    
    %    Loop to read the 562 fit images from the 44 with processing by 
    %    means of the Takeda method simulation, Cell 1.
    file = strcat('C:\Location on your Computer\NEW WORK\Interpoladas\Hong Wu\Image_Hong_', num2str(m), '.mat');
    
    %    Loop to read the 660 fit images from the 40 with processing by 
    %    means of the Takeda method simulation, Cell 2.
    %file = strcat('C:\Location on your Computer\NEW WORK\Interpoladas\Ullah\Image_Ullah_', num2str(m), '.mat');        

    EstrAux = load(file);
    fitimage = EstrAux.mat;
    %   Identification of images generated by 2DLegendre_Fit  
    aux_0 = fitimage;
    aux_1 = imresize(aux_0, [A B]);
    FitImage(:,:,m) = aux_1;
    %    Images with edges and zeros in NaN variables
    
    Aux = FitImage(:,:,m);
    for i=1:A
        for j=1:B
            if Aux(i,j) == 0
                NanImages(i,j,m) = NaN;
            else
                NanImages(i,j,m) = Aux(i,j);
            end
        end
    end
end

%    Routine to generate a contour plot for the 3D object with size MxMxM,
%    where each contour or level is one of the images generated by 
%    "2DLegendre_Fit". We bear in mind that this plot does not allow us to
%    observe the inner parts of te 3D object, only the surface and outer
%    edge of each level are observed.

figure;
hold on
for ii = 1 : M
    zz = (J) - (ii);
    [~,h] = contour(X,Y,NanImages(:,:,ii),30);  % plot contour at the bottom
    h.ContourZLevel = zz;
end
hold off
view(3);
xlabel('X [nm]'); ylabel('Y [nm]'); zlabel('Z [nm]', 'Rotation',0)

%    Instruction to visualize the 3D model (vision) of the SARS-CoV-2 viral
%    particle, which allows us to manipulate the object by rotating it with
%    respect to the XYZ axes independently, observe transversal cuts fixing
%    the XY-YZ-XZ planes along the axes X,Y and Z. Likewise, it allows us 
%    observing the inner parts with intensity filters, highlighting the 
%    interest zones with colors manually injected by the user. The latter
%    can be seen in the supplementary material in Videos 1 and 2 for Cell 1 
%    and Cell 2, respectively.

volumeViewer(NanImages);
 
