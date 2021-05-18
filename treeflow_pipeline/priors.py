import numpy as np
import treeflow
import treeflow.priors
import treeflow_pipeline.model


def get_prior_from_spec(dist_spec, probs):
    if dist_spec == "fixed":
        return dist_spec
    else:
        dist_name, quantiles = next(iter(dist_spec.items()))
        dist_class = treeflow_pipeline.model.dists[dist_name]
        params, res = treeflow.priors.get_params_for_quantiles(
            dist_class, quantiles, probs=probs
        )
        assert res.success
        return {dist_name: {key: value.item() for key, value in params.items()}}


def get_priors_from_spec(spec):
    probs = spec.pop("probs")
    clock_model, clock_spec = next(iter(spec.pop("clock").items()))
    if clock_model == "relaxed_lognormal_conjugate":
        (
            conjugate_prior_params,
            clock_res,
        ) = treeflow.priors.get_params_for_quantiles_lognormal_conjugate(
            clock_spec["cov"], clock_spec["mean"], probs=probs
        )
        assert clock_res["loc"].success
        assert clock_res["precision"].success
        clock_params = { f"rate_{param}": prior for param, prior in conjugate_prior_params.items() }
    else:
        raise ValueError(f"Clock model not implemented: {clock_model}")

    def get_model_part(model_part_spec):
        if model_part_spec == "none":
            return model_part_spec
        else:
            model_part, model_part_param_spec = next(iter(model_part_spec.items()))
            return {
                model_part: {
                    param_key: get_prior_from_spec(param_spec, probs)
                    for param_key, param_spec in model_part_param_spec.items()
                }
            }

    return {
        "clock": {
            clock_model: {
                param: {  # TODO: Refactor this with above function
                    dist_name: {
                        key: value.item() for key, value in prior_params.items()
                    }
                    for dist_name, prior_params in param_prior.items()
                }
                for param, param_prior in clock_params.items()
            }
        },
        **{
            model_part: get_model_part(model_part_spec)
            for model_part, model_part_spec in spec.items()
        },
    }


def get_free_params_from_model_spec(model_spec):
    clock_model = next(iter(model_spec["clock"]))
    if clock_model == "relaxed_lognormal_conjugate":
        clock_params = ["rate_loc", "rate_precision"]
    else:
        raise ValueError(f"Unknown clock model: {clock_model}")

    def get_free_params_from_model_spec_part(model_spec_part):
        if model_spec_part == "none":
            return []
        else:
            part_params = next(iter(model_spec_part.values()))
            return [
                param
                for param, param_model in part_params.items()
                if param_model != "fixed"
            ]

    model_spec_parts = [model_spec[part] for part in ["tree", "substitution", "site"]]
    other_params = [
        param
        for model_spec_part in model_spec_parts
        for param in get_free_params_from_model_spec_part(model_spec_part)
    ]
    return other_params + clock_params
