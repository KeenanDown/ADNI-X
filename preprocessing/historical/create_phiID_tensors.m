function create_phiID_tensors(BOLD_FILE, ELPH_LOCATION, PARCELLATION_LOCATION, PARCELLATION_NAME)
    try
        disp(BOLD_FILE)

        % This will save a .mat file at the BOLD_FILE location.

        % Calculate the number of parcels.
        NUMBER_OF_PARCELS = length(unique(niftiread(PARCELLATION_LOCATION))) - 1;

        addpath(ELPH_LOCATION);
        addpath(fullfile(ELPH_LOCATION, 'PhiID'));

        % We now parcellate the data.
        disp('Now performing the parcellation.')
        parcellated_scan = parcellate(BOLD_FILE, PARCELLATION_LOCATION, PARCELLATION_NAME);
        disp('Parcellation complete.')

        % Use file parts.
        file_parts = fileparts(BOLD_FILE);

        % Get the location to output just the parcellation output.
        output_location_parcellation = fullfile(file_parts, 'ParcellationMatrix.mat');

        % Save the file.
        save(output_location_parcellation, 'parcellated_scan');

        % Initialise the PhiIDTensor.
        PhiIDTensor = zeros(NUMBER_OF_PARCELS, NUMBER_OF_PARCELS, 16);

        % Compute the PhiID for every off-diagonal pair.
        for i = 1:NUMBER_OF_PARCELS
            disp(['Now computing for parcel ', num2str(i)])
            for j = 1:NUMBER_OF_PARCELS
                if i == j
                    continue   % skip diagonal entries
                end

                PhiIDTensorStruct = PhiIDFull(parcellated_scan([i, j], :));

                PhiIDTensor(i,j,1)  = PhiIDTensorStruct.rtr;
                PhiIDTensor(i,j,2)  = PhiIDTensorStruct.rtx;
                PhiIDTensor(i,j,3)  = PhiIDTensorStruct.rty;
                PhiIDTensor(i,j,4)  = PhiIDTensorStruct.rts;
                PhiIDTensor(i,j,5)  = PhiIDTensorStruct.xtr;
                PhiIDTensor(i,j,6)  = PhiIDTensorStruct.xtx;
                PhiIDTensor(i,j,7)  = PhiIDTensorStruct.xty;
                PhiIDTensor(i,j,8)  = PhiIDTensorStruct.xts;
                PhiIDTensor(i,j,9)  = PhiIDTensorStruct.ytr;
                PhiIDTensor(i,j,10) = PhiIDTensorStruct.ytx;
                PhiIDTensor(i,j,11) = PhiIDTensorStruct.yty;
                PhiIDTensor(i,j,12) = PhiIDTensorStruct.yts;
                PhiIDTensor(i,j,13) = PhiIDTensorStruct.str;
                PhiIDTensor(i,j,14) = PhiIDTensorStruct.stx;
                PhiIDTensor(i,j,15) = PhiIDTensorStruct.sty;
                PhiIDTensor(i,j,16) = PhiIDTensorStruct.sts;

                disp("Success.")
            end
        end

        disp("Finished assigning PhiID values.")

        % -----------------  (no symmetry-forcing step) -----------------

        % Get the location for an output.
        output_location = fullfile(file_parts, 'PhiIDTensor.mat');

        % Save the file.
        save(output_location, 'PhiIDTensor');
        disp(['Saved PhiIDTensor at ', output_location])

    catch ME
        disp(['Error: ', ME.message]);
    end
end
