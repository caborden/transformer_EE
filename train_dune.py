import json
import os

from transformer_ee.train import MVtrainer

print(os.getcwd())
with open("git/cborden_transformer_EE/transformer_EE/transformer_ee/config/mod_input_DUNE_atmo-E-th.json", encoding="UTF-8", mode="r") as f:
    input_d = json.load(f)

#input_d["data_path"]="/exp/dune/app/users/cborden/atmoDUNE/anatree_hd_AV_2dot5_random_sum_300k.root"
input_d["data_path"]= "/exp/dune/app/users/cborden/atmoDUNE/ana_tree_hd_999.root"
# input_d["model"]["name"] = "Transformer_EE_v4"
input_d["model"]["kwargs"]["nhead"] = 3
input_d["model"]["epochs"] = 5
input_d["model"]["kwargs"]["num_layers"] = 5
#input_d["optimizer"]["name"] = "sgd"
input_d["optimizer"]["name"] = "Adam"
input_d["optimizer"]["kwargs"]["lr"] = 0.001
#input_d["optimizer"]["kwargs"]["momentum"] = 0.9
#input_d["save_path"] = "/home/tthakore/save/model/DUNE_atmo/"
input_d["save_path"] = "/home/cborden/save/model/atmoDUNE/10_23_2024"
# input_d["weight"] = {"name": "FlatSpectraWeights", "kwargs": {"maxweight": 5, "minweight": 0.2}}

with open("git/cborden_transformer_EE/transformer_EE/transformer_ee/config/mod_input_DUNE_atmo-E-th.json", encoding="UTF-8", mode="w") as f:
    json.dump(input_d, f)


my_trainer = MVtrainer(input_d)
my_trainer.train()
my_trainer.eval()