import pytest
import torch
from transformers import AutoTokenizer, EsmModel
import os


RUN_ESM_MODEL_TESTS = os.getenv("RUN_ESM_MODEL_TESTS") == "1"

@pytest.mark.skipif(
    not RUN_ESM_MODEL_TESTS,
    reason="Set RUN_ESM_MODEL_TESTS=1 to run model-loading integration tests.",
)
def test_esm2_padding_behavior():
    """
    Verify if ESM-2 last_hidden_state contains garbage in padded positions.
    If it does, the manual masking in `esm_embedder.py` is NECESSARY.
    If it is already zeroed, the manual masking is redundant but harmless.
    """
    model_name = "facebook/esm2_t6_8M_UR50D"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = EsmModel.from_pretrained(model_name)
    model.eval()

    # Create a batch with one short and one long sequence
    seqs = ["MK", "M" * 20] 
    toks = tokenizer(seqs, return_tensors="pt", padding=True)
    
    with torch.no_grad():
        outputs = model(**toks)
        hs = outputs.last_hidden_state # [Batch, seq_len, hidden_dim]

    # Inspection
    # Sequence 1 is "MK" (length 2). 
    # Tokens: [BOS, M, K, EOS, PAD, PAD, ...]
    # We expect indices 4+ to be PAD.
    
    attention_mask = toks["attention_mask"][0] # 1 for non-pad
    
    # Check values at padded positions for the first sequence
    padded_vectors = hs[0, attention_mask == 0, :]
    
    # If means are non-zero, then manual masking IS required.
    non_zero_count = torch.count_nonzero(padded_vectors)
    total_elements = padded_vectors.numel()
    
    print(f"\nPadded elements: {total_elements}")
    print(f"Non-zero elements: {non_zero_count}")
    
    if non_zero_count > 0:
        print("CONCLUSION: ESM-2 DOES produce non-zero embeddings for PAD tokens.")
        print("The manual masking logic IS REQUIRED.")
    else:
        print("CONCLUSION: ESM-2 output is already zeroed for PAD tokens.")
        print("The manual masking logic is REDUNDANT.")
    
    assert hs.shape[0] == len(seqs)
    assert padded_vectors.numel() > 0
