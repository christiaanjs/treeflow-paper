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
        return {dist_name: params}


def get_priors_from_spec(spec):
    probs = spec.pop("probs")
    clock_model, clock_spec = next(iter(spec.pop("clock").items()))
    if clock_model == "relaxed_lognormal_conjugate":
        (
            clock_params,
            clock_res,
        ) = treeflow.priors.get_params_for_quantiles_lognormal_conjugate(
            clock_spec["cov"], clock_spec["mean"], probs=probs
        )
        assert clock_res["mean"].success
        assert clock_res["precision"].success
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
        "clock": {clock_model: clock_params},
        **{
            model_part: get_model_part(model_part_spec)
            for model_part, model_part_spec in spec.items()
        },
    }
