clc
clear all
close all

% Matlab script for outputing compression metrics for a given 16bit DICOM
% image and a target MSE value. Here Adaptive quantization with Diagonal
% Traversal is used. File path and target MSE value should be given as
% input parameters.


% Start the timer
tic

%input Parameters
fileName = ' ';
targetMSE = 225;


s = dir(fileName);
orgSize = s.bytes;

img = dicomread(fileName);
class = class(img);

if strcmp(class, 'uint16')

    dicomClass = @uint16;

elseif strcmp(class, 'int16')

    dicomClass = @int16;

else
    error('Error: Invalid DICOM class. Only supports uint16 and int16');
end


metadata = dicominfo(fileName);

w = 'db1';
n = 1;
[C,S] = wavedec2(img,n,w);

Csort = sort(abs(C(:)));

A1 = appcoef2(C,S,w,1);
[H1 V1 D1] = detcoef2('all',C,S,1);


iterationMat = [];
MSEMat = [];
peaksnrMat = [];
absDiffMat = [];
compressionRatioMat = [];
encMatSizeMat = [];


A1filt = [];
H1filt = [];
V1filt = [];
D1filt = [];
MSE = 0;
iteration = 0;
q = 0;
exitOuterLoop = false;
reEval = true;
startThresh = 0.20;
startQunt = 0;
EndQunt = 0;


if targetMSE<=20

    startQunt = 20;
    EndQunt = 1;

elseif targetMSE<=50

    startQunt = 30;
    EndQunt = 10;

elseif targetMSE<=100

    startQunt = 40;
    EndQunt = 20;

elseif targetMSE<=200

    startQunt = 50;
    EndQunt = 25;

else
    startThresh = 0.10;
    startQunt = 60;
    EndQunt = 20;

end

endThresh = 0.99;


for i = 1:4

    MSEMat = [];
    middleThresh = (startThresh + endThresh)/2;

    keep = [startThresh , middleThresh ];

    for i = 1:length(keep)
        keepValue = keep(i);

        thresh = Csort(round((1-keepValue)*length(Csort)));
        q = EndQunt;

        indA1 = abs(A1)>thresh;
        indH1 = abs(H1)>thresh;
        indV1 = abs(V1)>thresh;
        indD1 = abs(D1)>thresh;

        A1filt = A1.*indA1;
        H1filt = H1.*indH1;
        V1filt = V1.*indV1;
        D1filt = D1.*indD1;

        A1filt1 = round(A1filt/EndQunt);
        H1filt1 = round(H1filt/EndQunt);
        V1filt1 = round(V1filt/EndQunt);
        D1filt1 = round(D1filt/EndQunt);

        A1filtNew = A1filt1 * EndQunt;
        V1filtNew = V1filt1 * EndQunt;
        H1filtNew = H1filt1 * EndQunt;
        D1FiltNew = D1filt1 * EndQunt;

        reconstructedImage = dicomClass(idwt2(A1filtNew,H1filtNew,V1filtNew,D1FiltNew,w));
        MSEMat = [MSEMat immse(img,reconstructedImage)];

    end

    if MSEMat(1) >= targetMSE && MSEMat(2) <= targetMSE

        endThresh = middleThresh;

    else

        startThresh = middleThresh;

    end

end

startThresh = floor(startThresh*100) / 100;
endThresh = ceil(endThresh*100) / 100;


for n= startThresh: 0.001 : endThresh

    keep = [n];
    iteration = n;

    thresh = Csort(round((1-keep)*length(Csort)));

    q = startQunt;

    indA1 = abs(A1)>thresh;
    indH1 = abs(H1)>thresh;
    indV1 = abs(V1)>thresh;
    indD1 = abs(D1)>thresh;

    A1filt = A1.*indA1;
    H1filt = H1.*indH1;
    V1filt = V1.*indV1;
    D1filt = D1.*indD1;

    A1filt1 = round(A1filt/startQunt);
    H1filt1 = round(H1filt/startQunt);
    V1filt1 = round(V1filt/startQunt);
    D1filt1 = round(D1filt/startQunt);

    A1filtNew = A1filt1 * startQunt;
    V1filtNew = V1filt1 * startQunt;
    H1filtNew = H1filt1 * startQunt;
    D1FiltNew = D1filt1 * startQunt;


    reconstructedImage = dicomClass(idwt2(A1filtNew,H1filtNew,V1filtNew,D1FiltNew,w));

    MSE_atStartQ = immse(img,reconstructedImage);

    A1filt2 = round(A1filt/EndQunt);
    H1filt2 = round(H1filt/EndQunt);
    V1filt2 = round(V1filt/EndQunt);
    D1filt2 = round(D1filt/EndQunt);

    A1filtNew = A1filt2 * EndQunt;
    V1filtNew = V1filt2 * EndQunt;
    H1filtNew = H1filt2 * EndQunt;
    D1FiltNew = D1filt2 * EndQunt;

    reconstructedImage = dicomClass(idwt2(A1filtNew,H1filtNew,V1filtNew,D1FiltNew,w));
    MSE_atEndQ = immse(img,reconstructedImage);

    if MSE_atStartQ>=targetMSE && targetMSE>=MSE_atEndQ

        for p= startQunt: -1 : EndQunt

            A1filtQ = round(A1filt/q);
            H1filtQ = round(H1filt/q);
            V1filtQ = round(V1filt/q);
            D1filtQ = round(D1filt/q);

            A1filtNew = A1filtQ * q;
            V1filtNew = V1filtQ * q;
            H1filtNew = H1filtQ * q;
            D1FiltNew = D1filtQ * q;

            reconstructedImage = dicomClass(idwt2(A1filtNew,H1filtNew,V1filtNew,D1FiltNew,w));

            MSE = immse(img,reconstructedImage);

            if MSE>=targetMSE-0.1 && MSE<=targetMSE+0.1
                exitOuterLoop = true;
                reEval = false;
                break;

            end

            q = q - 1;

        end

    end

    if exitOuterLoop==true
        break;
    end

end

if reEval == true && EndQunt > 10

    startQunt = EndQunt;
    EndQunt = EndQunt - 10;

    for n= startThresh: 0.001 : endThresh

        keep = [n];
        iteration = n;

        thresh = Csort(round((1-keep)*length(Csort)));

        q = startQunt;

        indA1 = abs(A1)>thresh;
        indH1 = abs(H1)>thresh;
        indV1 = abs(V1)>thresh;
        indD1 = abs(D1)>thresh;

        A1filt = A1.*indA1;
        H1filt = H1.*indH1;
        V1filt = V1.*indV1;
        D1filt = D1.*indD1;

        A1filt1 = round(A1filt/startQunt);
        H1filt1 = round(H1filt/startQunt);
        V1filt1 = round(V1filt/startQunt);
        D1filt1 = round(D1filt/startQunt);

        A1filtNew = A1filt1 * startQunt;
        V1filtNew = V1filt1 * startQunt;
        H1filtNew = H1filt1 * startQunt;
        D1FiltNew = D1filt1 * startQunt;

        reconstructedImage = dicomClass(idwt2(A1filtNew,H1filtNew,V1filtNew,D1FiltNew,w));

        MSE_atStartQ = immse(img,reconstructedImage);

        A1filt2 = round(A1filt/EndQunt);
        H1filt2 = round(H1filt/EndQunt);
        V1filt2 = round(V1filt/EndQunt);
        D1filt2 = round(D1filt/EndQunt);

        A1filtNew = A1filt2 * EndQunt;
        V1filtNew = V1filt2 * EndQunt;
        H1filtNew = H1filt2 * EndQunt;
        D1FiltNew = D1filt2 * EndQunt;

        reconstructedImage = dicomClass(idwt2(A1filtNew,H1filtNew,V1filtNew,D1FiltNew,w));

        MSE_atEndQ = immse(img,reconstructedImage);

        if MSE_atStartQ>=targetMSE && targetMSE>=MSE_atEndQ

            for p= startQunt: -1 : EndQunt

                A1filtQ = round(A1filt/q);
                H1filtQ = round(H1filt/q);
                V1filtQ = round(V1filt/q);
                D1filtQ = round(D1filt/q);

                A1filtNew = A1filtQ * q;
                V1filtNew = V1filtQ * q;
                H1filtNew = H1filtQ * q;
                D1FiltNew = D1filtQ * q;

                reconstructedImage = dicomClass(idwt2(A1filtNew,H1filtNew,V1filtNew,D1FiltNew,w));

                MSE = immse(img,reconstructedImage);

                if MSE>=targetMSE-0.1 && MSE<=targetMSE+0.1
                    exitOuterLoop = true;
                    reEval = false;
                    break;

                end
                q = q - 1;

            end

        end

        if exitOuterLoop==true
            break;
        end

    end

end

A1filt = round(A1filt/q);
H1filt = round(H1filt/q);
V1filt = round(V1filt/q);
D1filt = round(D1filt/q);

reshapedA1 = reshape(A1filt, 1 , []);

VHD = [V1filt H1filt D1filt];

diagTraversedVHD = traverseDiagonal(VHD);

rleEncodedVHD = run_length_encode(diagTraversedVHD,'encode');

CombinedMatrix = [reshapedA1 rleEncodedVHD];

rowsA1 = size(A1filt,1);

coloumnsA1  = size(A1filt,2);

currentDir = pwd;

[~, fileName, ~] = fileparts(fileName);

CompFolder = strcat(fileName , '_Compressed');

if exist(CompFolder, 'dir')
    rmdir(CompFolder, 's');  % Remove existing folder and its contents
end

mkdir(CompFolder);

cd(CompFolder);

save("metadata.mat","q","metadata","rowsA1","coloumnsA1");

writematrix(CombinedMatrix);

gzip('CombinedMatrix.txt')

openfile = fopen("CombinedMatrix.txt.gz","r");

data = fread(openfile);

OpenDatafileID = fopen('ImgData1.bin', 'w');

fwrite(OpenDatafileID, data);

% Close the file
fclose(OpenDatafileID);
fclose(openfile);

delete('CombinedMatrix.txt',"CombinedMatrix.txt.gz");



%---------------Decompression----------------------

gunzip('ImgData1.bin')

CombinedMatrixRead = readmatrix('ImgData1');

load('metadata.mat');

colsComMat = size(CombinedMatrixRead,2);

A1ext = CombinedMatrix(1 , 1:(coloumnsA1*rowsA1));

VHDext = CombinedMatrix(1 , (coloumnsA1*rowsA1 + 1):colsComMat);

A1deQuantized = reshape(A1ext, rowsA1 , coloumnsA1) * q;

decodedVHD = run_length_encode(VHDext,'decode');

VHDTraversed = inverseTraverseDiagonal(decodedVHD, rowsA1 , coloumnsA1*3);

V1deQuantized = VHDTraversed(1:rowsA1 , 1:coloumnsA1) * q;

H1deQuantized = VHDTraversed(1:rowsA1 , (coloumnsA1 + 1) : coloumnsA1*2) * q;

D1deQuantized = VHDTraversed(1:rowsA1 , (coloumnsA1*2 + 1) : coloumnsA1*3) * q;

if isfield(metadata, 'BitsAllocated')
    bitsAllocated = metadata.BitsAllocated;
    if bitsAllocated == 16
        if isfield(metadata, 'PixelRepresentation') && metadata.PixelRepresentation == 1
            dicomClass = @int16;
        else
            dicomClass = @uint16;
        end
    else
        error('Error: Invalid DICOM class. Only supports uint16 and int16');
    end
else
    error('Error: Bit depth not found');
end

reconstructedImage = dicomClass(idwt2(A1deQuantized,H1deQuantized,V1deQuantized,D1deQuantized,w));

MSECHK = immse(img,reconstructedImage);

diffImage = imabsdiff(img,reconstructedImage);
meanDiff = mean(diffImage(:));

[peaksnr, snr] = psnr(reconstructedImage, img);

encMat = dir('ImgData1.bin');
encMatSize = encMat.bytes;

fileNameMetadata = 'metadata.mat';
metadatadir = dir(fileNameMetadata);
metdadataSize = metadatadir.bytes;

encMatSize = encMat.bytes + metdadataSize;

compressionRatio = orgSize/encMatSize;

MSECHK =  round(vpa(sym(MSECHK), 4),4);
meanDiff =  round(vpa(sym(meanDiff), 4),4);
peaksnr =  round(vpa(sym(peaksnr), 4),4);
compressionRatio = round(vpa(sym(compressionRatio), 4),4);

RMSE = sqrt(MSECHK);
SSIM = ssim(img,reconstructedImage);

RMSE =  round(vpa(sym(RMSE), 4),4);

targetMSE
MSECHK
RMSE
meanDiff
peaksnr
SSIM
compressionRatio

cd(currentDir)

% Stop the timer
elapsedTime = toc

%%
%RLE
function output = run_length_encode(data, type)
if strcmp(type, 'encode')
    output = [];
    i = 1;
    while i <= length(data)
        count = 1;
        while (i+count) <= length(data) && data(i+count) == data(i)
            count = count + 1;
        end
        output = [output, data(i), count];
        i = i + count;
    end
elseif strcmp(type, 'decode')
    output = [];
    for i = 1:2:length(data)-1
        element = data(i);
        count = data(i+1);
        output = [output, repmat(element, 1, count)];
    end
else
    error('Invalid encoding type');
end
end


%%
%traverseDiagonal

function output = traverseDiagonal(matrix)

[rows, cols] = size(matrix);

output = [];

for diag = 1 : rows+cols-1

    start_col = max([1, (diag - rows + 1)]);

    count = min( [diag, (cols - start_col + 1 ), rows] );

    line = [];

    for i = 1 : count

        element = matrix(min(rows, diag) - i  + 1, start_col + i - 1);

        line = [line , element];

    end

    output = [output , line];

end

end

%%
%inverseTraverseDiagonal

function output = inverseTraverseDiagonal(matrix, rows, cols)

rowPos = 1;
output = zeros(rows,cols);

for diag = 1 : rows+cols-1

    start_col = max([1, (diag - rows + 1)]);

    count = min( [diag, (cols - start_col + 1 ), rows] );

    for i = 1 : count

        output(min(rows, diag) - i  + 1, start_col + i - 1) = matrix(rowPos);

        rowPos = rowPos + 1;

    end

end

end
