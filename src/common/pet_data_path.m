function p = pet_data_path
% Tries to find path to MR raw data

if exist('SIRF_PET_DATA_PATH', 'var')
    p = SIRF_PET_DATA_PATH;
else
    SIRF_PATH = getenv('SIRF_PATH');
    if ~isempty(SIRF_PATH)
        p = [SIRF_PATH '/data/examples/PET'];
    else
        p = './';
    end
end
end