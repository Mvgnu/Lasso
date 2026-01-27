"""Interactive BGC Visualization for Pipeline Results.

Provides per-GBK visualization with:
- DNA sequence with CDS annotations
- Candidate positions with hover-to-highlight
- Selection capability for downstream processing
- Detailed sequence view with annotations
"""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any, Dict, List, Optional


def _escape(s: str) -> str:
    """HTML-escape a string."""
    return html.escape(str(s) if s else "")


def generate_pipeline_visualization_html(
    results_json: List[Dict[str, Any]],
    title: str = "Pipeline Results",
) -> str:
    """Generate interactive HTML visualization for pipeline results."""
    return generate_gradio_visualization_html(
        results_json,
        initial_selection=None,
        title=title,
        export_enabled=True,
    )


def generate_inline_visualization_html(
    results_json: List[Dict[str, Any]],
    title: str = "Pipeline Results",
) -> str:
    """Generate a simplified inline HTML visualization for Gradio embedding.
    
    This version uses an iframe with srcdoc to properly sandbox the visualization.
    """
    # Get the full HTML
    full_html = generate_pipeline_visualization_html(results_json, title)
    
    # Escape the HTML for use in srcdoc attribute
    escaped_html = full_html.replace("&", "&amp;").replace('"', "&quot;")
    
    # Create an iframe wrapper that will render the visualization
    return f'''<div style="width:100%;border-radius:8px;overflow:hidden;background:#f0f2f5;">
    <iframe 
        srcdoc="{escaped_html}"
        style="width:100%;height:800px;border:none;border-radius:8px;"
        sandbox="allow-scripts allow-same-origin allow-downloads"
    ></iframe>
</div>'''


def save_pipeline_visualization(
    results_json: List[Dict[str, Any]],
    output_path: Path,
    title: str = "Pipeline Results",
) -> Path:
    """Save visualization to HTML file."""
    html_content = generate_pipeline_visualization_html(results_json, title)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html_content)
    return output_path


def generate_gradio_visualization_html(
    results_json: List[Dict[str, Any]],
    initial_selection: Optional[List[str]] = None,
    title: str = "Pipeline Results",
    export_enabled: bool = False,
) -> str:
    """Generate HTML visualization optimized for Gradio embedding with selection sync.

    Features:
    - Clean slate/gray styling matching workbench theme
    - postMessage bridge for real-time selection sync with Gradio
    - Optional export buttons for standalone HTML
    - Pre-populated selection from initial_selection parameter

    Args:
        results_json: Pipeline results in JSON format
        initial_selection: List of candidate_ids that should be pre-selected
        title: Title for the visualization
        export_enabled: Show export buttons in the selection bar
    """
    data_json = json.dumps(results_json, indent=None)
    initial_sel_json = json.dumps(initial_selection or [])
    export_controls_html = ""
    if export_enabled:
        export_controls_html = """
            <div class="selection-actions">
                <button class="btn btn-primary btn-sm" onclick="exportSelected('fasta')">Export FASTA</button>
                <button class="btn btn-primary btn-sm" onclick="exportSelected('json')">Export JSON</button>
                <button class="btn btn-secondary btn-sm" onclick="clearSelection()">Clear</button>
            </div>
        """

    return f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{_escape(title)}</title>
    <style>
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #f5f5f5;
            color: #1f2937;
            line-height: 1.5;
            padding: 16px;
        }}

        /* Header - slate theme */
        .header {{
            background: #475569;
            color: white;
            padding: 16px 20px;
            border-radius: 8px;
            margin-bottom: 16px;
        }}
        .header h1 {{ font-size: 1.25rem; margin-bottom: 8px; font-weight: 600; }}
        .header-stats {{
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
            font-size: 0.85rem;
            opacity: 0.9;
        }}
        .header-stat {{ display: flex; align-items: center; gap: 5px; }}
        .header-stat strong {{ font-size: 1.1rem; }}

        /* Controls */
        .controls {{
            background: white;
            padding: 12px 16px;
            border-radius: 8px;
            margin-bottom: 16px;
            display: flex;
            gap: 16px;
            align-items: center;
            flex-wrap: wrap;
            border: 1px solid #e5e7eb;
        }}
        .control-group {{ display: flex; align-items: center; gap: 8px; }}
        .control-group label {{ font-size: 0.8rem; color: #64748b; font-weight: 500; }}
        select, input[type="range"] {{
            padding: 5px 8px;
            border: 1px solid #d1d5db;
            border-radius: 4px;
            font-size: 0.8rem;
            background: white;
        }}
        select:focus, input:focus {{ outline: 2px solid #64748b; outline-offset: 1px; }}
        .threshold-display {{
            font-weight: 600;
            color: #475569;
            min-width: 36px;
            font-size: 0.85rem;
        }}

        /* Selection info bar */
        .selection-info {{
            background: #f0fdf4;
            border: 1px solid #86efac;
            padding: 10px 16px;
            border-radius: 6px;
            margin-bottom: 16px;
            display: none;
            align-items: center;
            gap: 12px;
            font-size: 0.85rem;
        }}
        .selection-info.visible {{ display: flex; }}
        .selection-count {{ font-weight: 600; color: #166534; }}
        .selection-actions {{ display: flex; gap: 8px; margin-left: auto; }}
        .btn {{
            padding: 6px 10px;
            border: 1px solid transparent;
            border-radius: 4px;
            font-size: 0.75rem;
            cursor: pointer;
            font-weight: 600;
        }}
        .btn-primary {{ background: #0ea5a0; color: white; }}
        .btn-secondary {{ background: #e2e8f0; color: #1f2937; }}
        .btn-sm {{ padding: 4px 8px; }}

        /* Legend */
        .legend {{
            background: white;
            padding: 10px 14px;
            border-radius: 6px;
            margin-bottom: 16px;
            display: flex;
            gap: 14px;
            flex-wrap: wrap;
            align-items: center;
            font-size: 0.75rem;
            border: 1px solid #e5e7eb;
        }}
        .legend-title {{ font-weight: 600; color: #374151; }}
        .legend-item {{ display: flex; align-items: center; gap: 4px; color: #6b7280; }}
        .legend-color {{ width: 12px; height: 12px; border-radius: 2px; }}

        /* BGC Panel */
        .bgc-panel {{
            background: white;
            border-radius: 8px;
            margin-bottom: 12px;
            overflow: hidden;
            border: 1px solid #e5e7eb;
        }}
        .bgc-header {{
            background: #64748b;
            color: white;
            padding: 12px 16px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            cursor: pointer;
            font-size: 0.9rem;
        }}
        .bgc-header.lasso {{ background: #059669; }}
        .bgc-header:hover {{ filter: brightness(1.05); }}
        .bgc-title {{
            display: flex;
            align-items: center;
            gap: 8px;
            font-weight: 600;
        }}
        .bgc-title .source {{ opacity: 0.8; font-size: 0.8rem; font-weight: normal; }}
        .badge {{
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 0.65rem;
            font-weight: 600;
            text-transform: uppercase;
            background: rgba(255,255,255,0.2);
        }}
        .bgc-meta {{
            font-size: 0.75rem;
            opacity: 0.9;
            display: flex;
            gap: 12px;
        }}

        .bgc-body {{ padding: 0; }}
        .bgc-panel.collapsed .bgc-body {{ display: none; }}

        /* Segment Visualization */
        .segment-section {{
            padding: 14px 16px;
            border-bottom: 1px solid #e5e7eb;
        }}
        .segment-title {{
            font-size: 0.8rem;
            font-weight: 600;
            color: #64748b;
            margin-bottom: 10px;
        }}
        .segment-scroll {{
            overflow-x: auto;
            overflow-y: hidden;
            padding-bottom: 8px;
        }}
        .segment-canvas {{
            position: relative;
            height: 100px;
            background: #f8fafc;
            border: 1px solid #e5e7eb;
            border-radius: 4px;
            min-width: 100%;
            display: inline-block;
        }}

        /* DNA backbone */
        .dna-backbone {{
            position: absolute;
            left: 16px;
            right: 16px;
            top: 40px;
            height: 2px;
            background: #cbd5e1;
        }}
        .dna-backbone::before {{
            content: "5'";
            position: absolute;
            left: -14px;
            top: -5px;
            font-size: 9px;
            color: #94a3b8;
        }}
        .dna-backbone::after {{
            content: "3'";
            position: absolute;
            right: -14px;
            top: -5px;
            font-size: 9px;
            color: #94a3b8;
        }}

        /* CDS features */
        .cds-feature {{
            position: absolute;
            height: 20px;
            top: 31px;
            border-radius: 3px;
            cursor: pointer;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 7px;
            color: white;
            font-weight: 500;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            padding: 0 3px;
            transition: all 0.1s;
            opacity: 0.85;
        }}
        .cds-feature:hover {{ opacity: 1; transform: scaleY(1.1); z-index: 10; }}
        .cds-feature.highlighted {{ opacity: 1; transform: scaleY(1.15); z-index: 15; box-shadow: 0 0 6px currentColor; }}
        .cds-peptidase {{ background: #dc2626; }}
        .cds-cyclase {{ background: #0ea5a0; }}
        .cds-b_protein {{ background: #2563eb; }}
        .cds-precursor {{ background: #059669; }}
        .cds-transporter {{ background: #d97706; }}
        .cds-regulator {{ background: #0891b2; }}
        .cds-other {{ background: #6b7280; }}

        /* Candidate markers */
        .candidate-marker {{
            position: absolute;
            height: 14px;
            top: 60px;
            border-radius: 2px;
            cursor: pointer;
            transition: all 0.1s;
            background: #94a3b8;
            border: 1px solid #94a3b8;
        }}
        .candidate-marker:hover, .candidate-marker.highlighted {{
            transform: scaleY(1.3);
            z-index: 20;
        }}
        .candidate-marker.highlighted {{ box-shadow: 0 0 8px currentColor; }}
        .candidate-marker.known {{ border-color: #dc2626; }}
        .candidate-marker.novel {{ border-color: #475569; }}
        .candidate-marker.selected {{
            outline: 2px solid #475569;
            outline-offset: 1px;
        }}

        /* Ruler */
        .ruler {{
            position: absolute;
            bottom: 4px;
            left: 16px;
            right: 16px;
            height: 18px;
            font-size: 8px;
            color: #94a3b8;
        }}
        .ruler-tick {{ position: absolute; bottom: 10px; width: 1px; height: 5px; background: #cbd5e1; }}
        .ruler-label {{ position: absolute; bottom: 0; transform: translateX(-50%); white-space: nowrap; }}

        /* Two column layout */
        .content-columns {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 16px;
            padding: 16px;
        }}
        @media (max-width: 900px) {{ .content-columns {{ grid-template-columns: 1fr; }} }}

        /* CDS List */
        .annotation-section, .candidates-section {{
            background: #f8fafc;
            border-radius: 6px;
            overflow: hidden;
            border: 1px solid #e5e7eb;
        }}
        .section-header {{
            background: #f1f5f9;
            padding: 8px 12px;
            font-weight: 600;
            font-size: 0.8rem;
            color: #374151;
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid #e5e7eb;
        }}
        .section-header .count {{
            background: #475569;
            color: white;
            padding: 1px 6px;
            border-radius: 8px;
            font-size: 0.7rem;
        }}

        .annotation-list, .candidate-list {{
            max-height: 280px;
            overflow-y: auto;
        }}

        /* CDS row */
        .cds-row {{
            padding: 8px 12px;
            border-bottom: 1px solid #f1f5f9;
            cursor: pointer;
            transition: background 0.1s;
            font-size: 0.8rem;
        }}
        .cds-row:hover, .cds-row.highlighted {{ background: #e0f2fe; }}
        .cds-row:last-child {{ border-bottom: none; }}
        .cds-row-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 3px;
        }}
        .cds-type-badge {{
            padding: 1px 5px;
            border-radius: 2px;
            font-size: 0.65rem;
            font-weight: 600;
            color: white;
        }}
        .cds-product {{ color: #64748b; font-size: 0.75rem; }}
        .cds-location {{ color: #94a3b8; font-size: 0.7rem; }}

        /* Candidate row */
        .candidate-row {{
            padding: 10px 12px;
            border-bottom: 1px solid #f1f5f9;
            cursor: pointer;
            transition: background 0.1s;
        }}
        .candidate-row:hover, .candidate-row.highlighted {{ background: #faf5ff; }}
        .candidate-row.selected {{
            background: #f0fdf4;
            border-left: 3px solid #22c55e;
        }}
        .candidate-row:last-child {{ border-bottom: none; }}

        .candidate-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 5px;
        }}
        .candidate-header-left {{
            display: flex;
            align-items: center;
            gap: 6px;
        }}
        .candidate-checkbox {{
            width: 14px;
            height: 14px;
            cursor: pointer;
            accent-color: #22c55e;
        }}
        .candidate-id {{
            font-weight: 600;
            font-size: 0.75rem;
            color: #374151;
        }}
        .candidate-score {{
            font-weight: 700;
            font-size: 0.85rem;
        }}
        .score-high {{ color: #059669; }}
        .score-medium {{ color: #d97706; }}
        .score-low {{ color: #dc2626; }}

        .candidate-sequence {{
            font-family: 'SF Mono', Monaco, Consolas, monospace;
            font-size: 10px;
            background: white;
            padding: 5px 7px;
            border-radius: 3px;
            word-break: break-all;
            line-height: 1.4;
            border: 1px solid #e5e7eb;
            color: #374151;
        }}

        .score-bar {{
            height: 4px;
            background: #e5e7eb;
            border-radius: 3px;
            overflow: hidden;
            margin: 6px 0 4px;
        }}
        .score-fill {{
            height: 100%;
            border-radius: 3px;
        }}

        .frame-label {{
            position: absolute;
            left: 4px;
            font-size: 9px;
            color: #64748b;
        }}
        .frame-line {{
            position: absolute;
            left: 16px;
            right: 16px;
            height: 1px;
            background: #e2e8f0;
            opacity: 0.7;
        }}

        .candidate-meta {{
            display: flex;
            gap: 10px;
            margin-top: 5px;
            font-size: 0.7rem;
            color: #6b7280;
        }}

        .type-badge {{
            padding: 1px 5px;
            border-radius: 2px;
            font-size: 0.65rem;
            font-weight: 600;
        }}
        .type-badge.known {{ background: #fef2f2; color: #b91c1c; }}
        .type-badge.novel {{ background: #f0fdf4; color: #15803d; }}
        .type-badge.inverted {{ background: #ffedd5; color: #c2410c; }}

        /* Tooltip */
        .tooltip {{
            position: fixed;
            background: #1e293b;
            color: white;
            padding: 10px 14px;
            border-radius: 6px;
            font-size: 0.75rem;
            z-index: 1000;
            pointer-events: none;
            max-width: 360px;
            box-shadow: 0 4px 16px rgba(0,0,0,0.25);
            display: none;
        }}
        .tooltip.visible {{ display: block; }}
        .tooltip-title {{ font-weight: 600; margin-bottom: 5px; font-size: 0.8rem; }}
        .tooltip-sequence {{
            font-family: 'SF Mono', Monaco, Consolas, monospace;
            font-size: 9px;
            background: rgba(255,255,255,0.1);
            padding: 5px;
            border-radius: 3px;
            word-break: break-all;
            margin-top: 6px;
        }}
        .tooltip-row {{
            display: flex;
            justify-content: space-between;
            margin: 2px 0;
            font-size: 0.7rem;
        }}
        .tooltip-label {{ opacity: 0.7; }}

        /* Empty state */
        .empty-state {{
            padding: 24px;
            text-align: center;
            color: #94a3b8;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>{_escape(title)}</h1>
        <div class="header-stats" id="header-stats"></div>
    </div>

    <div class="controls">
        <div class="control-group">
            <label>Min Similarity:</label>
            <input type="range" id="threshold" min="0" max="100" value="0">
            <span class="threshold-display" id="threshold-display">0.00</span>
        </div>
        <div class="control-group">
            <label>Zoom:</label>
            <input type="range" id="zoom" min="0.5" max="500" step="0.1" value="1">
            <span class="threshold-display" id="zoom-display">1.0x</span>
        </div>
        <div class="control-group">
            <label>Show:</label>
            <select id="limit-select">
                <option value="10">Top 10</option>
                <option value="20">Top 20</option>
                <option value="50" selected>Top 50</option>
                <option value="100">Top 100</option>
                <option value="250">Top 250</option>
                <option value="500">Top 500</option>
                <option value="all">All</option>
            </select>
        </div>
        <div class="control-group">
            <label>Type:</label>
            <select id="type-filter">
                <option value="all">All</option>
                <option value="novel">Novel only</option>
                <option value="known">Known only</option>
            </select>
        </div>
        <div class="control-group">
            <label>Layout:</label>
            <select id="lane-mode">
                <option value="frame" selected>By frame</option>
                <option value="packed">Packed</option>
            </select>
        </div>
    </div>

    <div class="selection-info" id="selection-info">
        <span>Selected: <span class="selection-count" id="selection-count">0</span> candidates</span>
        {export_controls_html}
    </div>

    <div class="legend">
        <span class="legend-title">CDS:</span>
        <div class="legend-item"><div class="legend-color" style="background:#dc2626"></div>Peptidase</div>
        <div class="legend-item"><div class="legend-color" style="background:#0ea5a0"></div>Cyclase</div>
        <div class="legend-item"><div class="legend-color" style="background:#2563eb"></div>B-protein</div>
        <div class="legend-item"><div class="legend-color" style="background:#059669"></div>Precursor</div>
        <div class="legend-item"><div class="legend-color" style="background:#d97706"></div>Transporter</div>
        <div class="legend-item"><div class="legend-color" style="background:#0891b2"></div>Regulator</div>
        <div class="legend-item"><div class="legend-color" style="background:#6b7280"></div>Other</div>
        <span class="legend-title" style="margin-left:16px;">Score:</span>
        <div class="legend-item"><div class="legend-gradient"></div>Low → High</div>
        <span class="legend-title" style="margin-left:12px;">Border:</span>
        <div class="legend-item"><div class="legend-color" style="background:transparent;border:2px solid #dc2626;"></div>Known</div>
        <div class="legend-item"><div class="legend-color" style="background:transparent;border:2px solid #475569;"></div>Novel</div>
    </div>

    <div id="bgc-panels"></div>

    <div class="tooltip" id="tooltip"></div>

    <script>
const DATA = {data_json};
let selectedCandidates = new Set({initial_sel_json});
let currentThreshold = 0;
let currentLimit = 50;
let currentTypeFilter = 'all';
let currentZoom = 1;
let zoomRenderTimer = null;
let laneMode = 'frame';
const FRAME_ORDER = ['+0', '+1', '+2', '-1', '-2', '-3'];

// CDS classification
const CDS_KEYWORDS = {{
    peptidase: ['b2 protein', 'peptidase', 'protease'],
    cyclase: ['cyclase', 'isopeptide', 'transglutaminase'],
    b_protein: ['asparagine synthase', 'asparagine synthetase', 'b1 protein', 'b3 protein'],
    precursor: ['precursor', 'lasso'],
    transporter: ['transporter', 'abc', 'export'],
    regulator: ['regulator', 'transcription']
}};

function classifyCds(product) {{
    const p = (product || '').toLowerCase();
    for (const [type, keywords] of Object.entries(CDS_KEYWORDS)) {{
        for (const kw of keywords) {{
            if (p.includes(kw)) return type;
        }}
    }}
    return 'other';
}}

function escapeHtml(s) {{
    const d = document.createElement('div');
    d.textContent = s || '';
    return d.innerHTML;
}}

function formatNumber(n) {{
    return n.toLocaleString();
}}

function getScore(c) {{
    const score = c.combined_score || c.top_n_mean_similarity || c.best_similarity || 0;
    return isFinite(score) ? score : 0;
}}

function scoreToColor(score) {{
    const clamped = Math.max(0, Math.min(1, score));
    const hue = 120 * clamped;
    return `hsl(${{hue}}, 70%, 45%)`;
}}

function frameLabel(c) {{
    if (c.frame === null || c.frame === undefined) {{
        return 'unk';
    }}
    const frame = parseInt(c.frame, 10);
    if (Number.isNaN(frame)) {{
        return 'unk';
    }}
    return frame >= 0 ? `+${{frame}}` : `${{frame}}`;
}}

function getFrameKeys(candidates) {{
    const seen = new Set();
    (candidates || []).forEach(c => {{
        const label = frameLabel(c);
        if (label !== 'unk') {{
            seen.add(label);
        }}
    }});
    const ordered = FRAME_ORDER.filter(f => seen.has(f));
    if (ordered.length === 0 && seen.size > 0) {{
        return Array.from(seen).sort();
    }}
    return ordered.length ? ordered : ['unk'];
}}

// Notify parent (Gradio) of selection changes
function pushSelectionToReceiver() {{
    try {{
        const doc = (window.top && window.top.document) ? window.top.document : window.parent?.document;
        if (!doc) return;
        let receiver = doc.getElementById('viz-selection-receiver');
        if (receiver && receiver.tagName !== 'TEXTAREA' && receiver.tagName !== 'INPUT') {{
            receiver = receiver.querySelector('textarea, input');
        }}
        if (!receiver) {{
            receiver = doc.querySelector('#viz-selection-receiver textarea, #viz-selection-receiver input');
        }}
        if (!receiver) return;
        receiver.value = JSON.stringify(Array.from(selectedCandidates));
        receiver.dispatchEvent(new Event('input', {{ bubbles: true }}));
        receiver.dispatchEvent(new Event('change', {{ bubbles: true }}));
    }} catch (e) {{
        // Best-effort; postMessage is the primary path.
    }}
}}

function notifySelectionChange() {{
    pushSelectionToReceiver();
    const selectionArray = Array.from(selectedCandidates);
    try {{
        window.parent.postMessage({{
            type: 'candidate_selection',
            candidates: selectionArray
        }}, '*');
        if (window.top && window.top !== window.parent) {{
            window.top.postMessage({{
                type: 'candidate_selection',
                candidates: selectionArray
            }}, '*');
        }}
    }} catch (e) {{
        console.log('Could not post message to parent:', e);
    }}
    updateSelectionUI();
}}

// Initialize
function init() {{
    updateStats();
    renderAllPanels();
    setupEventListeners();
    updateSelectionUI();
    // Notify parent of initial selection
    if (selectedCandidates.size > 0) {{
        notifySelectionChange();
    }}
}}

function updateStats() {{
    let totalBgcs = DATA.length;
    let lassoBgcs = DATA.filter(r => r.is_lasso).length;
    let totalCandidates = DATA.reduce((sum, r) => sum + (r.candidates?.length || 0), 0);
    let visibleCandidates = DATA.reduce((sum, r) => sum + getFilteredCandidates(r.candidates).length, 0);

    document.getElementById('header-stats').innerHTML = `
        <div class="header-stat"><strong>${{totalBgcs}}</strong> BGCs</div>
        <div class="header-stat"><strong>${{lassoBgcs}}</strong> Lasso</div>
        <div class="header-stat"><strong>${{formatNumber(totalCandidates)}}</strong> Total ORFs</div>
        <div class="header-stat"><strong>${{formatNumber(visibleCandidates)}}</strong> Visible</div>
    `;
}}

function getFilteredCandidates(candidates) {{
    if (!candidates) return [];

    let filtered = candidates.filter(c => {{
        const score = c.combined_score || c.top_n_mean_similarity || c.best_similarity || 0;
        if (score < currentThreshold) return false;

        const isKnown = Boolean(c.is_known_precursor);
        if (currentTypeFilter === 'novel' && isKnown) return false;
        if (currentTypeFilter === 'known' && !isKnown) return false;

        return true;
    }});

    // Sort by score
    filtered.sort((a, b) => {{
        const scoreA = a.combined_score || a.top_n_mean_similarity || a.best_similarity || 0;
        const scoreB = b.combined_score || b.top_n_mean_similarity || b.best_similarity || 0;
        return scoreB - scoreA;
    }});

    // Apply limit
    if (currentLimit !== 'all' && filtered.length > currentLimit) {{
        filtered = filtered.slice(0, currentLimit);
    }}

    return filtered;
}}

let scrollState = {{}};

function captureScrollState() {{
    const container = document.getElementById('bgc-panels');
    if (!container) return;
    scrollState = {{}};
    const nodes = container.querySelectorAll('[data-scroll-key]');
    nodes.forEach((node) => {{
        const key = node.getAttribute('data-scroll-key');
        if (!key) return;
        const max = node.scrollWidth - node.clientWidth;
        const leftRatio = max > 0 ? (node.scrollLeft / max) : 0;
        scrollState[key] = {{
            leftRatio,
            scrollTop: node.scrollTop,
        }};
    }});
}}

function restoreScrollState() {{
    const container = document.getElementById('bgc-panels');
    if (!container) return;
    const nodes = container.querySelectorAll('[data-scroll-key]');
    nodes.forEach((node) => {{
        const key = node.getAttribute('data-scroll-key');
        const saved = key ? scrollState[key] : null;
        if (!saved) return;
        const max = node.scrollWidth - node.clientWidth;
        node.scrollLeft = max > 0 ? (saved.leftRatio * max) : 0;
        node.scrollTop = saved.scrollTop;
    }});
}}

function renderAllPanels() {{
    const container = document.getElementById('bgc-panels');
    const scrollY = window.scrollY;
    const scrollX = window.scrollX;
    captureScrollState();
    container.innerHTML = '';

    if (DATA.length === 0) {{
        container.innerHTML = `<div class="empty-state">No BGC data to display</div>`;
        return;
    }}

    DATA.forEach((bgc, idx) => {{
        container.appendChild(renderBgcPanel(bgc, idx));
    }});
    applySelectionStyles();
    requestAnimationFrame(() => {{
        window.scrollTo(scrollX, scrollY);
        restoreScrollState();
        requestAnimationFrame(restoreScrollState);
        setTimeout(restoreScrollState, 0);
    }});
}}

function renderBgcPanel(bgc, idx) {{
    const panel = document.createElement('div');
    panel.className = 'bgc-panel';
    panel.id = `bgc-${{idx}}`;

    const isLasso = bgc.is_lasso;
    const cds = bgc.annotated_cds || [];
    const candidates = getFilteredCandidates(bgc.candidates);
    const segmentLength = bgc.segment_length_nt || 1;
    const sourceFile = bgc.source_file ? bgc.source_file.split('/').pop() : '';

    panel.innerHTML = `
        <div class="bgc-header ${{isLasso ? 'lasso' : ''}}" onclick="togglePanel(${{idx}})">
            <div class="bgc-title">
                <span>${{escapeHtml(bgc.record_id)}}</span>
                ${{sourceFile ? `<span class="source">(${{escapeHtml(sourceFile)}})</span>` : ''}}
                <span class="badge">${{isLasso ? 'LASSO' : 'OTHER'}}</span>
                <span class="badge">${{candidates.length}} candidates</span>
            </div>
            <div class="bgc-meta">
                <span>${{formatNumber(segmentLength)}} bp</span>
                <span>${{cds.length}} CDS</span>
            </div>
        </div>
        <div class="bgc-body">
            <div class="segment-section">
                <div class="segment-title">Genomic Segment</div>
                <div class="segment-scroll" data-scroll-key="segment-${{idx}}">
                    ${{renderSegmentCanvas(bgc, candidates, idx)}}
                </div>
            </div>
            <div class="content-columns">
                <div class="annotation-section">
                    <div class="section-header">
                        <span>CDS Annotations</span>
                        <span class="count">${{cds.length}}</span>
                    </div>
                    <div class="annotation-list" data-scroll-key="cds-${{idx}}">
                        ${{renderCdsList(cds, idx)}}
                    </div>
                </div>
                <div class="candidates-section">
                    <div class="section-header">
                        <span>Candidates</span>
                        <span class="count">${{candidates.length}}</span>
                    </div>
                    <div class="candidate-list" data-scroll-key="candidates-${{idx}}">
                        ${{renderCandidateList(candidates, idx)}}
                    </div>
                </div>
            </div>
        </div>
    `;

    return panel;
}}

function renderSegmentCanvas(bgc, candidates, bgcIdx) {{
    const length = bgc.segment_length_nt || 1;
    const cds = bgc.annotated_cds || [];

    const minWidth = 700;
    const maxWidth = 1800;
    const idealWidth = Math.max(minWidth, Math.min(maxWidth, length * 0.04));
    const scale = (idealWidth / length) * currentZoom;
    const canvasWidth = length * scale + 32;
    const baseHeight = 110;
    const laneHeight = 18;
    const markerBaseTop = 68;
    const packedOverlap = 4;
    const markerHeight = laneMode === 'packed' ? 12 : 14;
    const laneStep = laneMode === 'packed' ? (markerHeight - packedOverlap) : laneHeight;
    const frameKeys = laneMode === 'frame' ? getFrameKeys(candidates) : [];
    const frameIndex = new Map();
    frameKeys.forEach((key, idx) => {{
        frameIndex.set(key, idx);
    }});
    const sublaneStep = 4;
    const maxSublane = 2;
    const sublaneMap = new Map();
    if (laneMode === 'frame') {{
        const frameGroups = {{}};
        candidates.forEach(c => {{
            const key = frameLabel(c);
            if (!frameGroups[key]) frameGroups[key] = [];
            frameGroups[key].push(c);
        }});
        Object.entries(frameGroups).forEach(([key, group]) => {{
            const laneEndsLocal = [];
            [...group]
                .sort((a, b) => (a.genomic_start || 0) - (b.genomic_start || 0))
                .forEach(c => {{
                    const start = Math.min(c.genomic_start || 0, c.genomic_end || 0);
                    const end = Math.max(c.genomic_start || 0, c.genomic_end || 0);
                    let lane = 0;
                    while (lane < laneEndsLocal.length && start < laneEndsLocal[lane]) {{
                        lane += 1;
                    }}
                    if (lane === laneEndsLocal.length) {{
                        laneEndsLocal.push(end);
                    }} else {{
                        laneEndsLocal[lane] = Math.max(laneEndsLocal[lane], end);
                    }}
                    sublaneMap.set(c.candidate_id, Math.min(lane, maxSublane));
                }});
        }});
    }}

    const laneMap = new Map();
    const laneEnds = [];
    const sortedCandidates = [...candidates].sort(
        (a, b) => (a.genomic_start || 0) - (b.genomic_start || 0)
    );
    sortedCandidates.forEach(c => {{
        const start = Math.min(c.genomic_start || 0, c.genomic_end || 0);
        const end = Math.max(c.genomic_start || 0, c.genomic_end || 0);
        let lane = 0;
        while (lane < laneEnds.length && start < laneEnds[lane]) {{
            lane += 1;
        }}
        if (lane === laneEnds.length) {{
            laneEnds.push(end);
        }} else {{
            laneEnds[lane] = Math.max(laneEnds[lane], end);
        }}
        laneMap.set(c.candidate_id, lane);
    }});
    const laneCount = laneMode === 'frame'
        ? Math.max(1, frameKeys.length)
        : Math.max(1, laneEnds.length);
    const canvasHeight = baseHeight + (laneCount - 1) * laneStep;

    let cdsHtml = '';
    cds.forEach((c, cdsIdx) => {{
        const type = classifyCds(c.product);
        const left = 16 + (c.start || 0) * scale;
        const width = Math.max(3, ((c.end || 0) - (c.start || 0)) * scale);
        const label = width > 25 ? (c.product || '').substring(0, Math.floor(width / 4)) : '';

        cdsHtml += `
            <div class="cds-feature cds-${{type}}"
                 id="cds-${{bgcIdx}}-${{cdsIdx}}"
                 style="left:${{left}}px; width:${{width}}px;"
                 data-cds='${{JSON.stringify(c).replace(/'/g, "&#39;")}}'
                 onmouseenter="highlightCds(${{bgcIdx}}, ${{cdsIdx}}); showTooltip(event, 'cds', this)"
                 onmouseleave="unhighlightCds(${{bgcIdx}}, ${{cdsIdx}}); hideTooltip()">
                ${{escapeHtml(label)}}
            </div>
        `;
    }});

    let candidateHtml = '';
    candidates.forEach((c, candIdx) => {{
        const isKnown = Boolean(c.is_known_precursor);
        const isSelected = selectedCandidates.has(c.candidate_id);
        const score = getScore(c);
        const color = scoreToColor(score);
        const left = 16 + (c.genomic_start || 0) * scale;
        const width = Math.max(2, ((c.genomic_end || 0) - (c.genomic_start || 0)) * scale);
        const lane = laneMode === 'frame'
            ? (frameIndex.get(frameLabel(c)) ?? 0)
            : (laneMap.get(c.candidate_id) || 0);
        const sublaneOffset = laneMode === 'frame'
            ? (sublaneMap.get(c.candidate_id) || 0) * sublaneStep
            : 0;
        const top = markerBaseTop + lane * laneStep + sublaneOffset;

        candidateHtml += `
            <div class="candidate-marker ${{isKnown ? 'known' : 'novel'}} ${{isSelected ? 'selected' : ''}}"
                 id="marker-${{bgcIdx}}-${{c.candidate_id}}"
                 data-candidate-id="${{escapeHtml(c.candidate_id)}}"
                 style="left:${{left}}px; width:${{width}}px; top:${{top}}px; height:${{markerHeight}}px; background:${{color}};"
                 data-candidate='${{JSON.stringify(c).replace(/'/g, "&#39;")}}'
                 onclick="selectCandidate(${{bgcIdx}}, '${{c.candidate_id}}')"
                 onmouseenter="highlightCandidate(${{bgcIdx}}, '${{c.candidate_id}}'); showTooltip(event, 'candidate', this)"
                 onmouseleave="unhighlightCandidate(${{bgcIdx}}, '${{c.candidate_id}}'); hideTooltip()">
            </div>
        `;
    }});

    let frameGuideHtml = '';
    if (laneMode === 'frame' && frameKeys.length) {{
        frameKeys.forEach((key, idx) => {{
            const y = markerBaseTop + idx * laneStep + 2;
            frameGuideHtml += `
                <div class="frame-label" style="top:${{y}}px;">${{key}}</div>
                <div class="frame-line" style="top:${{y + 6}}px;"></div>
            `;
        }});
    }}

    let rulerHtml = '';
    const numTicks = 5;
    for (let i = 0; i <= numTicks; i++) {{
        const bp = Math.round((length / numTicks) * i);
        const x = 16 + bp * scale;
        rulerHtml += `
            <div class="ruler-tick" style="left:${{x}}px;"></div>
            <div class="ruler-label" style="left:${{x}}px;">${{formatNumber(bp)}}</div>
        `;
    }}

    return `
        <div class="segment-canvas" style="width:${{canvasWidth}}px; height:${{canvasHeight}}px;">
            <div class="dna-backbone" style="left:16px; right:16px;"></div>
            ${{cdsHtml}}
            ${{frameGuideHtml}}
            ${{candidateHtml}}
            <div class="ruler">
                ${{rulerHtml}}
            </div>
        </div>
    `;
}}

function renderCdsList(cds, bgcIdx) {{
    if (!cds || cds.length === 0) {{
        return '<div class="empty-state">No CDS annotations</div>';
    }}

    return cds.map((c, idx) => {{
        const type = classifyCds(c.product);
        const typeColors = {{
            peptidase: '#dc2626', cyclase: '#0ea5a0', b_protein: '#2563eb',
            precursor: '#059669', transporter: '#d97706', regulator: '#0891b2', other: '#6b7280'
        }};

        return `
            <div class="cds-row" id="cds-row-${{bgcIdx}}-${{idx}}"
                 onmouseenter="highlightCdsFeature(${{bgcIdx}}, ${{idx}})"
                 onmouseleave="unhighlightCdsFeature(${{bgcIdx}}, ${{idx}})">
                <div class="cds-row-header">
                    <span class="cds-type-badge" style="background:${{typeColors[type]}}">${{type.toUpperCase().replace('_', ' ')}}</span>
                    <span class="cds-location">${{formatNumber(c.start || 0)}} - ${{formatNumber(c.end || 0)}}</span>
                </div>
                <div class="cds-product">${{escapeHtml(c.product || 'Unknown product')}}</div>
            </div>
        `;
    }}).join('');
}}

function renderCandidateList(candidates, bgcIdx) {{
    if (!candidates || candidates.length === 0) {{
        return '<div class="empty-state">No candidates match filters</div>';
    }}

    return candidates.map((c, idx) => {{
        const isKnown = Boolean(c.is_known_precursor);
        const isSelected = selectedCandidates.has(c.candidate_id);
        const isInverted = c.rule_orientation === 'core_leader';
        const score = getScore(c);
        const scoreColor = scoreToColor(score);
        const scorePct = Math.round(Math.max(0, Math.min(1, score)) * 100);
        const ruleScore = c.rule_score_raw;
        const frameTag = frameLabel(c);

        return `
            <div class="candidate-row ${{isSelected ? 'selected' : ''}}"
                 id="cand-row-${{bgcIdx}}-${{c.candidate_id}}"
                 data-candidate-id="${{escapeHtml(c.candidate_id)}}"
                 onclick="selectCandidateRow(event, ${{bgcIdx}}, '${{c.candidate_id}}')"
                 onmouseenter="highlightCandidateMarker(${{bgcIdx}}, '${{c.candidate_id}}')"
                 onmouseleave="unhighlightCandidateMarker(${{bgcIdx}}, '${{c.candidate_id}}')">
                <div class="candidate-header">
                    <div class="candidate-header-left">
                        <input type="checkbox" class="candidate-checkbox"
                               ${{isSelected ? 'checked' : ''}}
                               onclick="event.stopPropagation(); toggleCandidateSelection('${{c.candidate_id}}')">
                        <span class="candidate-id">${{escapeHtml(c.candidate_id.substring(0, 20))}}</span>
                        <span class="type-badge ${{isKnown ? 'known' : 'novel'}}">${{isKnown ? 'KNOWN' : 'NOVEL'}}</span>
                        ${{isInverted ? '<span class="type-badge inverted">CORE+LEADER</span>' : ''}}
                    </div>
                    <span class="candidate-score" style="color:${{scoreColor}}">${{score.toFixed(4)}}</span>
                    ${{ruleScore !== null && ruleScore !== undefined ? `<span class="candidate-score" style="color:#64748b;font-weight:600;">R ${{ruleScore.toFixed(2)}}</span>` : ''}}
                </div>
                <div class="score-bar"><div class="score-fill" style="width:${{scorePct}}%; background:${{scoreColor}};"></div></div>
                <div class="candidate-sequence">${{escapeHtml(c.protein_sequence)}}</div>
                <div class="candidate-meta">
                    <span>${{c.aa_length || '?'}} aa</span>
                    <span>${{formatNumber(c.genomic_start || 0)}} - ${{formatNumber(c.genomic_end || 0)}}</span>
                    <span>${{c.strand === 1 ? '+' : '-'}} strand</span>
                    ${{frameTag !== 'unk' ? `<span>frame ${{frameTag}}</span>` : ''}}
                    ${{c.best_match_id ? `<span>${{escapeHtml(c.best_match_id)}}</span>` : ''}}
                </div>
            </div>
        `;
    }}).join('');
}}

// Highlighting functions
function highlightCds(bgcIdx, cdsIdx) {{
    document.getElementById(`cds-${{bgcIdx}}-${{cdsIdx}}`)?.classList.add('highlighted');
    document.getElementById(`cds-row-${{bgcIdx}}-${{cdsIdx}}`)?.classList.add('highlighted');
}}

function unhighlightCds(bgcIdx, cdsIdx) {{
    document.getElementById(`cds-${{bgcIdx}}-${{cdsIdx}}`)?.classList.remove('highlighted');
    document.getElementById(`cds-row-${{bgcIdx}}-${{cdsIdx}}`)?.classList.remove('highlighted');
}}

function highlightCdsFeature(bgcIdx, cdsIdx) {{
    document.getElementById(`cds-${{bgcIdx}}-${{cdsIdx}}`)?.classList.add('highlighted');
}}

function unhighlightCdsFeature(bgcIdx, cdsIdx) {{
    document.getElementById(`cds-${{bgcIdx}}-${{cdsIdx}}`)?.classList.remove('highlighted');
}}

function highlightCandidate(bgcIdx, candId) {{
    document.getElementById(`marker-${{bgcIdx}}-${{candId}}`)?.classList.add('highlighted');
    document.getElementById(`cand-row-${{bgcIdx}}-${{candId}}`)?.classList.add('highlighted');
}}

function unhighlightCandidate(bgcIdx, candId) {{
    document.getElementById(`marker-${{bgcIdx}}-${{candId}}`)?.classList.remove('highlighted');
    document.getElementById(`cand-row-${{bgcIdx}}-${{candId}}`)?.classList.remove('highlighted');
}}

function highlightCandidateMarker(bgcIdx, candId) {{
    document.getElementById(`marker-${{bgcIdx}}-${{candId}}`)?.classList.add('highlighted');
}}

function unhighlightCandidateMarker(bgcIdx, candId) {{
    document.getElementById(`marker-${{bgcIdx}}-${{candId}}`)?.classList.remove('highlighted');
}}

function scrollToCandidateRow(bgcIdx, candId) {{
    const row = document.getElementById(`cand-row-${{bgcIdx}}-${{candId}}`);
    if (row) {{
        row.scrollIntoView({{ block: 'center', behavior: 'smooth' }});
    }}
}}

function scrollToCandidateMarker(bgcIdx, candId) {{
    const marker = document.getElementById(`marker-${{bgcIdx}}-${{candId}}`);
    if (marker) {{
        marker.scrollIntoView({{ block: 'center', behavior: 'smooth' }});
    }}
}}

function selectCandidate(bgcIdx, candId) {{
    toggleCandidateSelection(candId);
    scrollToCandidateRow(bgcIdx, candId);
}}

function selectCandidateRow(evt, bgcIdx, candId) {{
    if (evt && evt.target && evt.target.classList.contains('candidate-checkbox')) {{
        return;
    }}
    toggleCandidateSelection(candId);
    scrollToCandidateMarker(bgcIdx, candId);
}}

// Selection with notification
function toggleCandidateSelection(candId) {{
    if (selectedCandidates.has(candId)) {{
        selectedCandidates.delete(candId);
    }} else {{
        selectedCandidates.add(candId);
    }}
    notifySelectionChange();
    applySelectionStyles();
}}

function clearSelection() {{
    selectedCandidates.clear();
    notifySelectionChange();
    applySelectionStyles();
}}

function exportSelected(format) {{
    const selected = [];
    DATA.forEach(bgc => {{
        (bgc.candidates || []).forEach(c => {{
            if (selectedCandidates.has(c.candidate_id)) {{
                selected.push({{
                    ...c,
                    bgc_id: bgc.record_id,
                    source_file: bgc.source_file
                }});
            }}
        }});
    }});

    if (selected.length === 0) {{
        alert('No candidates selected');
        return;
    }}

    let content = '';
    let filename = '';

    if (format === 'fasta') {{
        content = selected.map(c =>
            `>${{c.candidate_id}} bgc=${{c.bgc_id}} score=${{getScore(c).toFixed(4)}}\\n${{c.protein_sequence}}`
        ).join('\\n');
        filename = 'selected_candidates.fasta';
    }} else {{
        content = JSON.stringify(selected, null, 2);
        filename = 'selected_candidates.json';
    }}

    const blob = new Blob([content], {{ type: 'text/plain' }});
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
}}

function applySelectionStyles() {{
    const rows = document.querySelectorAll('.candidate-row');
    rows.forEach(row => {{
        const candId = row.dataset.candidateId;
        if (!candId) return;
        const isSelected = selectedCandidates.has(candId);
        row.classList.toggle('selected', isSelected);
        const checkbox = row.querySelector('.candidate-checkbox');
        if (checkbox) checkbox.checked = isSelected;
    }});
    const markers = document.querySelectorAll('.candidate-marker');
    markers.forEach(marker => {{
        const candId = marker.dataset.candidateId;
        if (!candId) return;
        const isSelected = selectedCandidates.has(candId);
        marker.classList.toggle('selected', isSelected);
    }});
}}

function updateSelectionUI() {{
    const infoBar = document.getElementById('selection-info');
    const count = document.getElementById('selection-count');

    if (selectedCandidates.size > 0) {{
        infoBar.classList.add('visible');
        count.textContent = selectedCandidates.size;
    }} else {{
        infoBar.classList.remove('visible');
    }}
}}

// Receive selection updates from parent (Gradio)
window.addEventListener('message', (event) => {{
    if (event.data && event.data.type === 'set_selection') {{
        selectedCandidates = new Set(event.data.candidates || []);
        updateSelectionUI();
        applySelectionStyles();
    }}
}});

// Panel toggle
function togglePanel(idx) {{
    document.getElementById(`bgc-${{idx}}`)?.classList.toggle('collapsed');
}}

// Tooltip
function showTooltip(event, type, element) {{
    const tooltip = document.getElementById('tooltip');
    let html = '';

    if (type === 'cds') {{
        const data = JSON.parse(element.dataset.cds);
        const cdsType = classifyCds(data.product);
        html = `
            <div class="tooltip-title">${{escapeHtml(data.product || 'Unknown')}}</div>
            <div class="tooltip-row"><span class="tooltip-label">Type:</span> ${{cdsType.replace('_', ' ')}}</div>
            <div class="tooltip-row"><span class="tooltip-label">Location:</span> ${{formatNumber(data.start)}} - ${{formatNumber(data.end)}}</div>
            <div class="tooltip-row"><span class="tooltip-label">Strand:</span> ${{data.strand === 1 ? '+' : '-'}}</div>
            ${{data.locus_tag ? `<div class="tooltip-row"><span class="tooltip-label">Locus:</span> ${{escapeHtml(data.locus_tag)}}</div>` : ''}}
            ${{data.translation ? `<div class="tooltip-sequence">${{escapeHtml(data.translation.substring(0, 50))}}${{data.translation.length > 50 ? '...' : ''}}</div>` : ''}}
        `;
    }} else if (type === 'candidate') {{
        const data = JSON.parse(element.dataset.candidate);
        const isKnown = Boolean(data.is_known_precursor);
        html = `
            <div class="tooltip-title">${{escapeHtml(data.candidate_id)}}</div>
            <div class="tooltip-row"><span class="tooltip-label">Type:</span> ${{isKnown ? 'Known precursor' : 'Novel candidate'}}</div>
            <div class="tooltip-row"><span class="tooltip-label">Score:</span> ${{(data.combined_score || data.top_n_mean_similarity || data.best_similarity || 0).toFixed(4)}}</div>
            <div class="tooltip-row"><span class="tooltip-label">Length:</span> ${{data.aa_length}} aa</div>
            <div class="tooltip-row"><span class="tooltip-label">Location:</span> ${{formatNumber(data.genomic_start)}} - ${{formatNumber(data.genomic_end)}}</div>
            ${{data.best_match_id ? `<div class="tooltip-row"><span class="tooltip-label">Best match:</span> ${{escapeHtml(data.best_match_id)}}</div>` : ''}}
            <div class="tooltip-sequence">${{escapeHtml(data.protein_sequence)}}</div>
        `;
    }}

    tooltip.innerHTML = html;
    tooltip.classList.add('visible');

    const x = event.clientX + 12;
    const y = event.clientY + 12;
    tooltip.style.left = x + 'px';
    tooltip.style.top = y + 'px';

    const rect = tooltip.getBoundingClientRect();
    if (rect.right > window.innerWidth) {{
        tooltip.style.left = (x - rect.width - 24) + 'px';
    }}
    if (rect.bottom > window.innerHeight) {{
        tooltip.style.top = (y - rect.height - 24) + 'px';
    }}
}}

function hideTooltip() {{
    document.getElementById('tooltip').classList.remove('visible');
}}

// Event listeners
function setupEventListeners() {{
    const threshold = document.getElementById('threshold');
    const thresholdDisplay = document.getElementById('threshold-display');
    threshold.addEventListener('input', () => {{
        currentThreshold = threshold.value / 100;
        thresholdDisplay.textContent = currentThreshold.toFixed(2);
    }});
    threshold.addEventListener('change', () => {{
        updateStats();
        renderAllPanels();
    }});

    const zoom = document.getElementById('zoom');
    const zoomDisplay = document.getElementById('zoom-display');
    zoom.addEventListener('input', () => {{
        currentZoom = parseFloat(zoom.value);
        zoomDisplay.textContent = `${{currentZoom.toFixed(1)}}x`;
        if (zoomRenderTimer) {{
            clearTimeout(zoomRenderTimer);
        }}
        zoomRenderTimer = setTimeout(() => {{
            renderAllPanels();
            zoomRenderTimer = null;
        }}, 120);
    }});

    document.getElementById('limit-select').addEventListener('change', (e) => {{
        currentLimit = e.target.value === 'all' ? 'all' : parseInt(e.target.value);
        updateStats();
        renderAllPanels();
    }});

    document.getElementById('type-filter').addEventListener('change', (e) => {{
        currentTypeFilter = e.target.value;
        updateStats();
        renderAllPanels();
    }});

    document.getElementById('lane-mode').addEventListener('change', (e) => {{
        laneMode = e.target.value;
        renderAllPanels();
    }});
}}

// Initialize on load
window.addEventListener('DOMContentLoaded', init);
    </script>
</body>
</html>'''


def generate_gradio_iframe_html(
    results_json: List[Dict[str, Any]],
    initial_selection: Optional[List[str]] = None,
    title: str = "Pipeline Results",
    height: int = 700,
) -> str:
    """Generate an iframe wrapper for embedding the visualization in Gradio.

    This properly sandboxes the visualization while allowing postMessage communication.
    """
    inner_html = generate_gradio_visualization_html(
        results_json,
        initial_selection=initial_selection,
        title=title,
        export_enabled=False,
    )
    escaped_html = inner_html.replace("&", "&amp;").replace('"', "&quot;")

    return f'''<div style="width:100%;border-radius:8px;overflow:hidden;background:#f5f5f5;">
    <iframe
        id="viz-iframe"
        srcdoc="{escaped_html}"
        style="width:100%;height:{height}px;border:none;border-radius:8px;"
        sandbox="allow-scripts allow-same-origin allow-downloads"
    ></iframe>
</div>'''
