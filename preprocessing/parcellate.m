function ROI_activity = parcellate(BOLD_FILE, PARCELLATION_LOCATION, PARCELLATION_NAME)

    % Handle the case where the parcellation_name is not specified.
    if nargin < 3
        PARCELLATION_NAME = 'Parcellation'
    end

    % Specify the subject run to be parcellated.
    file_to_parcellate = BOLD_FILE;
    % Specify the location of the parcellation. NB this needs to be the
    % same dimensions as the volumes of the run.
    location_of_parcellation = PARCELLATION_LOCATION;

    % Import the subject data.
    disp("Importing the BOLD data for parcellation.")
    subject_data = niftiread(file_to_parcellate);
    % Import the parcellation data.
    disp("Importing the parcellation data.")
    parcellation_data = niftiread(location_of_parcellation);
    % Compute the number of parcels, removing the 0 parcels.
    number_of_parcels = length(unique(parcellation_data)) - 1;
    % Find the total number of volumes.
    disp("Finding the number of volumes.")
    number_of_volumes = niftiinfo(file_to_parcellate).ImageSize(4);

    % With thanks to Pedro for the following snippet.

    % Initialise the ROI activity matrix.
    ROI_activity = zeros([number_of_parcels, number_of_volumes]);
    % Iterating over the number of volumes,
    disp("Taking the mean across parcels.")
    for t = 1:number_of_volumes
        % Select the volume
        subject_at_t = subject_data(:,:,:,t);
        % For each parcel,
        for parcel = 1:number_of_parcels
            % Average the activity in that parcel and save.
            ROI_activity(parcel, t) = mean(subject_at_t(parcellation_data == parcel));
        end
    end
    disp("Finished.")

    % The following snippet from ChatGPT and perplexity.

    % Get the path of the BOLD_FILE
    [dirr, name, ext] = fileparts(BOLD_FILE);

    % Create the output file name
    % output_filename = char(strcat(dirr, filesep, name, '_', PARCELLATION_NAME, '_ROI_activity.mat'));
    % Unfortunately this seems to not notice the .nii before .nii.gz. So I need to add this.
    % ChatGPT gave me this:

    % Get the directory and name without extension
    [dirr, name, ext] = fileparts(BOLD_FILE);

    % Check if there is a double extension and remove both
    if strcmp(ext, '.gz')
        [~, name, ~] = fileparts(name);
    end

    % Create the output file name using the directory, base name without extensions, and PARCELLATION_NAME
    output_filename = char(strcat(dirr, filesep, name, '_', PARCELLATION_NAME, '_ROI_activity.mat'));


    % Save the ROI_activity variable as a MAT file in the same folder as the BOLD_FILE

    disp(output_filename)
    save(output_filename, 'ROI_activity');

    return
end
