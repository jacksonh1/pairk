import torch


class ESM_Model:
    """
    This is taken from the kibby conservation method and trivially modifed: DOI: 10.1093/bib/bbac599.
    see https://github.com/esbgkannan/kibby

    Class that loads a specified ESM model. Provides a method for encoding protein sequences.

    available models:
    - esm1b_t33_650M_UR50S
    - esm2_t6_8M_UR50D
    - esm2_t12_35M_UR50D
    - esm2_t30_150M_UR50D
    - esm2_t33_650M_UR50D (default)
    - esm2_t36_3B_UR50D


    Attributes
    ----------
    model_name : str
        the name of the model that was loaded.

    """

    def __init__(self, model_name: str = "esm2_t33_650M_UR50D"):
        self._load(model_name)

    def _load(self, model_name):
        import esm

        self.model_name = model_name
        self.model, alphabet = eval(f"esm.pretrained.{self.model_name}()")
        self.batch_converter = alphabet.get_batch_converter()
        self.model.eval()
        self.embed_dim = self.model._modules["layers"][0].embed_dim
        self.layers = sum(1 for i in self.model._modules["layers"])

    def encode(self, sequence, device="cuda", threads=1):
        """encode a protein sequence using the loaded model.

        Parameters
        ----------
        sequence : str
            the amino acid sequence to encode.
        device : str, optional
            whether to use a GPU via "cuda", or "cpu", by default "cuda"
        threads : int, optional
            The number of threads to use in pytorch, by default 1

        Returns
        -------
        torch.Tensor
            sequence embedding tensor
        """

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
                print("failed, trying CPU")
                return self.encode(sequence, device="cpu", threads=threads)
            else:
                return
