function data = precompute_MADE(model_ref)
% Pre-compute cobra model to elf model conversion.
%
% Author: Daniel Machado, 2013 

    data.elf_model = make_elf_model(model_ref);
    data.bounds_ref.lb = model_ref.lb;
    data.bounds_ref.ub = model_ref.ub;
end

