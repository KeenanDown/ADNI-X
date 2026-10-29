function create_phiID_tensors_v3(BOLD_FILE, ELPH_LOCATION, PARCELLATION_LOCATION, PARCELLATION_NAME)
%CREATE_PHIID_TENSORS_V3  Compute a full PhiID tensor for a parcellated fMRI scan
%                         using the **MMI** redundancy measure.
%
%   This version is identical to the original save-and-side-car workflow,
%   but calls PhiIDFull(..., 'MMI') so that the Minimum Mutual Information
%   (MMI) redundancy is used instead of the default CCS.
%
%   File-naming follows the BIDS-style pattern:
%       <BOLD-stem>_parc-<N>_desc-phiidtensor.mat / .json
%
%   Parameters
%   ----------
%   BOLD_FILE             : path to a 4-D pre-processed BOLD NIfTI
%   ELPH_LOCATION         : folder containing the PhiID toolbox
%   PARCELLATION_LOCATION : path to an integer-label NIfTI parcellation map
%   PARCELLATION_NAME     : short label for the parcellation atlas
%

    %----------------------------------------------------------------------
    % 0.  Settings
    %----------------------------------------------------------------------
    REDUNDANCY_MEASURE = 'MMI';   % <-- the only change required
    try
        disp(BOLD_FILE)

        %------------------------------------------------------------------
        % 1.  Determine number of parcels
        %------------------------------------------------------------------
        NUMBER_OF_PARCELS = length(unique(niftiread(PARCELLATION_LOCATION))) - 1; % subtract label 0

        %------------------------------------------------------------------
        % 2.  Add PhiID toolbox to the path
        %------------------------------------------------------------------
        addpath(ELPH_LOCATION);
        addpath(fullfile(ELPH_LOCATION, 'PhiID'));

        %------------------------------------------------------------------
        % 3.  Parcellate the data
        %------------------------------------------------------------------
        disp('Now performing the parcellation.');
        parcellated_scan = parcellate(BOLD_FILE, PARCELLATION_LOCATION, PARCELLATION_NAME);
        disp('Parcellation complete.');

        % Optionally save the parcellated matrix
        [scan_dir, scan_name, scan_ext] = fileparts(BOLD_FILE);
        if strcmpi(scan_ext, '.gz'), [~, scan_name] = fileparts(scan_name); end
        parcellation_out = fullfile(scan_dir, sprintf('%s_parcellationMatrix.mat', scan_name));
        save(parcellation_out, 'parcellated_scan');

        %------------------------------------------------------------------
        % 4.  Initialise PhiID tensor & compute all pairwise atoms
        %------------------------------------------------------------------
        PhiIDTensor = zeros(NUMBER_OF_PARCELS, NUMBER_OF_PARCELS, 16);

        for i = 1:NUMBER_OF_PARCELS
            fprintf('Now computing for parcel %d\n', i);
            for j = 1:NUMBER_OF_PARCELS
                if i == j, continue; end  % skip diagonal

                % --- key line: force MMI redundancy measure ----------------
                PhiIDTensorStruct = PhiIDFull( ...
                    parcellated_scan([i, j], :), [], REDUNDANCY_MEASURE);

                PhiIDTensor(i,j, 1) = PhiIDTensorStruct.rtr;
                PhiIDTensor(i,j, 2) = PhiIDTensorStruct.rtx;
                PhiIDTensor(i,j, 3) = PhiIDTensorStruct.rty;
                PhiIDTensor(i,j, 4) = PhiIDTensorStruct.rts;
                PhiIDTensor(i,j, 5) = PhiIDTensorStruct.xtr;
                PhiIDTensor(i,j, 6) = PhiIDTensorStruct.xtx;
                PhiIDTensor(i,j, 7) = PhiIDTensorStruct.xty;
                PhiIDTensor(i,j, 8) = PhiIDTensorStruct.xts;
                PhiIDTensor(i,j, 9) = PhiIDTensorStruct.ytr;
                PhiIDTensor(i,j,10) = PhiIDTensorStruct.ytx;
                PhiIDTensor(i,j,11) = PhiIDTensorStruct.yty;
                PhiIDTensor(i,j,12) = PhiIDTensorStruct.yts;
                PhiIDTensor(i,j,13) = PhiIDTensorStruct.str;
                PhiIDTensor(i,j,14) = PhiIDTensorStruct.stx;
                PhiIDTensor(i,j,15) = PhiIDTensorStruct.sty;
                PhiIDTensor(i,j,16) = PhiIDTensorStruct.sts;
            end
        end
        disp('Finished assigning PhiID values.');

        %------------------------------------------------------------------
        % 5.  Construct output filenames (BIDS-like)
        %------------------------------------------------------------------
        [~, stem, ext] = fileparts(BOLD_FILE);
        if strcmpi(ext, '.gz'), [~, stem] = fileparts(stem); end
        outStem = sprintf('%s_parc-%d_desc-phiidtensor', stem, NUMBER_OF_PARCELS);
        mat_out  = fullfile(scan_dir, [outStem, '.mat']);
        json_out = fullfile(scan_dir, [outStem, '.json']);

        %------------------------------------------------------------------
        % 6.  Save tensor and minimal provenance
        %------------------------------------------------------------------
        save(mat_out, 'PhiIDTensor');
        fprintf('Saved PhiIDTensor at %s\n', mat_out);

        meta.description        = 'Tensor of PhiID atoms';
        meta.created            = datestr(now, 'yyyy-mm-dd HH:MM:SS');
        meta.source_bold        = BOLD_FILE;
        meta.number_of_parcels  = NUMBER_OF_PARCELS;
        meta.redundancy_measure = REDUNDANCY_MEASURE;  % <-- extra field
        meta.fields_order       = {'rtr','rtx','rty','rts','xtr','xtx','xty','xts', ...
                                   'ytr','ytx','yty','yts','str','stx','sty','sts'};
        fid = fopen(json_out, 'w');
        fprintf(fid, '%s', jsonencode(meta, 'PrettyPrint', true));
        fclose(fid);
        fprintf('Wrote provenance side-car at %s\n', json_out);

    catch ME
        fprintf('Error: %s\n', ME.message);
    end
end
