from treeflow_pipeline.util import yaml_input, yaml_output
from treeflow_pipeline.model import generate_models_grid

grid_dict = yaml_input("config/data-models-grid.yaml")
data_models_dict = generate_models_grid(grid_dict)
yaml_output(data_models_dict, "config/data-models.yaml")
