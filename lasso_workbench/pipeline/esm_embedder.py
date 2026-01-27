"""
ESM-2 embedding wrapper
License: MIT
Author: Magnus Ohle
"""
import os
from typing import List, Optional

import torch
from transformers import AutoTokenizer, EsmModel


class ESM2Embedder:
    def __init__(
        self,
        model_name: Optional[str] = None,
        device: Optional[str] = None,
    ):
        if not model_name:
            raise ValueError("model_name must be provided")

        if device:
            self.device = device
        elif torch.cuda.is_available():
            self.device = "cuda"
        elif getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
            self.device = "mps"
        else:
            self.device = "cpu"

        self.tokenizer = AutoTokenizer.from_pretrained(
            model_name,
            do_lower_case=False,
        )
        # We only use last_hidden_state; disabling pooling avoids irrelevant pooler warnings.
        self.model = EsmModel.from_pretrained(
            model_name,
            add_pooling_layer=False,
        )
        self.model.eval().to(self.device)

    def embed_batch(self, sequences: List[str], batch_size: int = 1):
        """Embed sequences with padded batching (fast)."""
        import numpy as np

        seqs_all = list(sequences)
        bs = max(1, int(batch_size))
        outputs: list[np.ndarray] = []

        #https://docs.bioembeddings.com/v0.2.3/_modules/bio_embeddings/embed/esm_embedder.html ESM1 implementation used as reference
        with torch.inference_mode():
            for start in range(0, len(seqs_all), bs):
                seqs = seqs_all[start : start + bs]
                # tokenize batch of sequences in parallel, one sequence = one row in the batch
                toks = self.tokenizer(
                    seqs,
                    return_tensors="pt",
                    padding=True,
                    add_special_tokens=True,
                    return_special_tokens_mask=True,
                )
                # mask the special tokens to exclude padding/beginning of sequence/end of sequence tokens from the embedding
                # Get the value for special tokens based on the special_tokens_mask (1 for BOS/EOS/PAD, 0 for residues), remove key from dict
                special = toks.pop("special_tokens_mask")  
                # Get the value based on the attention mask is 1 for non‑pad tokens (including BOS/EOS)
                # 0 for pad tokens.
                attn = toks["attention_mask"]  
                # Construct the final mask (1 for residues, 0 for specials) 
                # [Batch, Seq_Len -> attn]*[Batch, Seq_Len -> 1-special]
                mask = attn * (1 - special)  
                # move inputs to device (gpu, mps, cpu)
                toks = {k: v.to(self.device) for k, v in toks.items()} 
                mask = mask.to(self.device)  # [Batch -> which sequence in the batch, Seq_Position -> which position in the sequence] now contains 1 for residues, 0 for specials
                # Get embedding per token (= amino acid) [Batch, Seq_Len, Hidden_Dim]
                hs = self.model(**toks).last_hidden_state  
                # Use mask, converted to (Batch, Seq_Len, 1 -> added for vector multiplication)
                # to multiply against all hidden states (individual token/amino acid embeddings)
                # BOS/EOS/PAD tokens get masked out (multiplied by 0)
                masked = hs * mask.unsqueeze(-1)  
                # Count valid residues per sequence (min 1)
                denom = mask.sum(dim=1, keepdim=True).clamp(min=1) 
                # Mean pooling: sum valid vectors / count of valid vectors 
                # calculate embedding per protein based on individual amino acid vectors for comparison 
                # (also: drop gradients, move to CPU, convert to NumPy, and force float32)
                pooled = (masked.sum(dim=1) / denom).detach().cpu().numpy().astype(np.float32, copy=False) 
                outputs.append(pooled)

        return np.vstack(outputs) # concatenate all batches


__all__ = ["ESM2Embedder"]
