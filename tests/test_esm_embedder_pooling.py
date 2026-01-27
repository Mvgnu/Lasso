import numpy as np
import torch

from lasso_workbench.pipeline.esm_embedder import ESM2Embedder


class _FakeModel:
    def __call__(self, **toks):
        input_ids = toks["input_ids"]
        bsz, seq_len = input_ids.shape
        hidden = 4
        hs = torch.zeros((bsz, seq_len, hidden), dtype=torch.float32)
        for i in range(bsz):
            for j in range(seq_len):
                hs[i, j, :] = float(j)
        return type("Out", (), {"last_hidden_state": hs})()

    def eval(self):
        return self

    def to(self, _device):
        return self


class _FakeTokenizer:
    def __call__(
        self,
        seqs,
        return_tensors="pt",
        padding=True,
        add_special_tokens=True,
        return_special_tokens_mask=True,
    ):
        max_len = max(len(s) for s in seqs) + 2  # BOS + EOS
        bsz = len(seqs)
        input_ids = torch.zeros((bsz, max_len), dtype=torch.long)
        attention_mask = torch.zeros((bsz, max_len), dtype=torch.long)
        special_mask = torch.zeros((bsz, max_len), dtype=torch.long)

        for i, seq in enumerate(seqs):
            length = len(seq) + 2
            attention_mask[i, :length] = 1
            special_mask[i, 0] = 1
            special_mask[i, length - 1] = 1
            if length < max_len:
                special_mask[i, length:] = 1

        return {
            "input_ids": input_ids,
            "attention_mask": attention_mask,
            "special_tokens_mask": special_mask,
        }


def test_pooling_excludes_special_tokens():
    embedder = ESM2Embedder.__new__(ESM2Embedder)
    embedder.device = "cpu"
    embedder.tokenizer = _FakeTokenizer()
    embedder.model = _FakeModel()

    pooled = embedder.embed_batch(["ACD"], batch_size=1)
    # Residues are positions 1..len(seq) => 1..3
    expected = np.array([[2.0, 2.0, 2.0, 2.0]], dtype=np.float32)
    np.testing.assert_allclose(pooled, expected)
