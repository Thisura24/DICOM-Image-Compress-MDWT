clc
clear all
close all

% Matlab script for outputing compression metrics for a given 16bit DICOM
% image and for threshold rate values from 0.2 to 0.99 with a step of 0.03.
% Here Adaptive quantization with Diagonal Traversal is used.
% File paths should be given as input parameters.

% Start the timer
tic

fileNames = {
    ' ', ' ', ' ', ' '
    };


for i = 1:numel(fileNames)

    % Clear workspace
    clearvars -except fileNames i;

    fileName = fileNames{i}

    s = dir(fileName);
    orgSize = s.bytes

    img = dicomread(fileName);
    class = class(img)

    metadata = dicominfo(fileName);

    if strcmp(class, 'uint16')

        dicomClass = @uint16;

    elseif strcmp(class, 'int16')

        dicomClass = @int16;

    else
        error('Error: Invalid DICOM class. Only supports uint16 and int16');
    end

    w = 'db1'
    n = 1;
    [C,S] = wavedec2(img,n,w);

    Csort = sort(abs(C(:)));

    A1 = appcoef2(C,S,w,1);
    [H1 V1 D1] = detcoef2('all',C,S,1);

    iterationMat = [];
    RMSEMat = [];
    peaksnrMat = [];
    absDiffMat = [];
    compressionRatioMat = [];
    ssimMat = [];
    MSEMargin = 0.05

    q = 0;

    for n= 0.2: 0.03 : 0.99

        keep1 = [n];

        thresh = Csort(round((1-keep1)*length(Csort)));

        indA1 = abs(A1)>thresh;
        indH1 = abs(H1)>thresh;
        indV1 = abs(V1)>thresh;
        indD1 = abs(D1)>thresh;

        A1filt = A1.*indA1;
        H1filt = H1.*indH1;
        V1filt = V1.*indV1;
        D1filt = D1.*indD1;

        for p= 5: 1 : 40

            q = p;

            indA1 = abs(A1)>thresh;
            indH1 = abs(H1)>thresh;
            indV1 = abs(V1)>thresh;
            indD1 = abs(D1)>thresh;

            A1filt = A1.*indA1;
            H1filt = H1.*indH1;
            V1filt = V1.*indV1;
            D1filt = D1.*indD1;

            reconstructedImage = dicomClass(idwt2(A1filt,H1filt,V1filt,D1filt,w));
            MSE1 = immse(img,reconstructedImage);

            A1filt = round(A1filt/q);
            H1filt = round(H1filt/q);
            V1filt = round(V1filt/q);
            D1filt = round(D1filt/q);

            A1filtNew = A1filt * q;
            V1reshaped = V1filt * q;
            H1reshaped = H1filt * q;
            D1reshaped = D1filt * q;

            reconstructedImage = dicomClass(idwt2(A1filtNew,H1reshaped,V1reshaped,D1reshaped,w));
            MSE2 = immse(img,reconstructedImage);

            if MSE2>= (MSE1 + MSEMargin*MSE1*n)

                break;

            end

        end

        reshapedA1 = reshape(A1filt, 1 , []);
        VHD = [V1filt H1filt D1filt];
        diagTraversedVHD = traverseDiagonal(VHD);
        rleEncodedVHD = run_length_encode(diagTraversedVHD,'encode');
        CombinedMatrix = [reshapedA1 rleEncodedVHD];

        writematrix(CombinedMatrix);

        rowsA1 = size(A1filt,1);

        coloumnsA1  = size(A1filt,2);

        writematrix(CombinedMatrix);

        save("data.mat","q","metadata","rowsA1","coloumnsA1");
        gzip('CombinedMatrix.txt');
        openfile = fopen("CombinedMatrix.txt.gz","r");
        data = fread(openfile);
        OpenDatafileID = fopen('ImgData1.bin', 'w');
        fwrite(OpenDatafileID, data);

        % Close the file
        fclose(OpenDatafileID);
        fclose(openfile);
        delete('CombinedMatrix.txt',"CombinedMatrix.txt.gz");
        rowsA1 = size(A1filt,1);
        coloumnsA1  = size(A1filt,2);
        colsComMat = size(CombinedMatrix,2);
        A1ext = CombinedMatrix(1 , 1:(coloumnsA1*rowsA1));
        VHDext = CombinedMatrix(1 , (coloumnsA1*rowsA1 + 1):colsComMat);
        A1filtNew = reshape(A1ext, rowsA1 , coloumnsA1) * q;
        decodedVHD = run_length_encode(VHDext,'decode');
        VHDTraversed = inverseTraverseDiagonal(decodedVHD, rowsA1 , coloumnsA1*3);

        V1filtNew = VHDTraversed(1:rowsA1 , 1:coloumnsA1) * q;

        H1filtNew = VHDTraversed(1:rowsA1 , (coloumnsA1 + 1) : coloumnsA1*2) * q;

        D1filtNew = VHDTraversed(1:rowsA1 , (coloumnsA1*2 + 1) : coloumnsA1*3) * q;

        reconstructedImage = dicomClass(idwt2(A1filtNew,H1filtNew,V1filtNew,D1filtNew,w));

        MSE = immse(img,reconstructedImage);

        diffImage = imabsdiff(img,reconstructedImage);
        meanDiff = mean(diffImage(:));
        [peaksnr, snr] = psnr(reconstructedImage, img);

        encMat = dir('ImgData1.bin');
        encMatSize = encMat.bytes;

        fileNameMetadata = 'data.mat';
        metadatadir = dir(fileNameMetadata);
        metdadataSize = metadatadir.bytes;

        encMatSize = encMat.bytes + metdadataSize;

        compressionRatio = orgSize/encMatSize;

        itr =  round(vpa(sym(n), 4),4);
        SSIM = ssim(img,reconstructedImage);
        RMSE = sqrt(MSE);
        RMSE =  round(vpa(sym(RMSE), 4),4);
        meanDiff =  round(vpa(sym(meanDiff), 4),4);
        peaksnr =  round(vpa(sym(peaksnr), 4),4);
        compressionRatio = round(vpa(sym(compressionRatio), 4),4);

        iterationMat = [iterationMat; itr];
        RMSEMat = [RMSEMat; RMSE];
        absDiffMat = [absDiffMat; meanDiff];
        peaksnrMat = [peaksnrMat; peaksnr];
        compressionRatioMat = [compressionRatioMat; compressionRatio];
        ssimMat = [ssimMat; SSIM];

    end

    iterationMat_numeric = double(iterationMat);
    RMSEMat_numeric = double(RMSEMat);
    absDiffMat_numeric = double(absDiffMat);
    peaksnrMat_numeric = double(peaksnrMat);
    ssimMat_numeric = double(ssimMat);
    compressionRatioMat_numeric = double(compressionRatioMat);

    % Write matrices to Excel file
    writematrix(iterationMat_numeric, 'data.xlsx', 'Sheet', i, 'Range', 'A1');
    writematrix(RMSEMat_numeric, 'data.xlsx', 'Sheet', i, 'Range', 'B1');
    writematrix(absDiffMat_numeric, 'data.xlsx', 'Sheet', i, 'Range', 'C1');
    writematrix(peaksnrMat_numeric, 'data.xlsx', 'Sheet', i, 'Range', 'D1');
    writematrix(ssimMat_numeric, 'data.xlsx', 'Sheet', i, 'Range', 'E1');
    writematrix(compressionRatioMat_numeric, 'data.xlsx', 'Sheet', i, 'Range', 'F1');

end

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

