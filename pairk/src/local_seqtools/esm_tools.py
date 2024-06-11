import torch

class ESM_Model:
    # esm1b_t33_650M_UR50S
    # esm2_t6_8M_UR50D
    # esm2_t12_35M_UR50D
    # esm2_t30_150M_UR50D
    # esm2_t33_650M_UR50D
    # esm2_t36_3B_UR50D

    def __init__(self, model_name: str = "esm2_t33_650M_UR50D"):
        self.load(model_name)

    def load(self, model_name):
        import esm

        self.model_name = model_name
        self.model, alphabet = eval(f"esm.pretrained.{self.model_name}()")
        self.batch_converter = alphabet.get_batch_converter()
        self.model.eval()
        self.embed_dim = self.model._modules["layers"][0].embed_dim
        self.layers = sum(1 for i in self.model._modules["layers"])

    def encode(self, sequence, device="cuda", threads=1):
        try:
            torch.cuda.empty_cache()
            torch.set_num_threads(threads)

            batch_labels, batch_strs, batch_tokens = self.batch_converter(
                [["", sequence]]
            )
            batch_tokens = batch_tokens.to(device)
            self.model = self.model.to(device)
            with torch.no_grad():
                results = self.model(
                    batch_tokens, repr_layers=[self.layers], return_contacts=False
                )
                results = results["representations"][self.layers].to("cpu")[0]
            return results
        except:
            if device != "cpu":
                return self.encode(sequence, device="cpu", threads=threads)
            else:
                return