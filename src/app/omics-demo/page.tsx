"use client";

import { useState } from "react";
import Link from "next/link";

// Sample RNA-seq data
const sampleData = [
  { gene: "TP53", lfc: 3.2, padj: 0.0001, baseMean: 5432 },
  { gene: "MYC", lfc: 2.8, padj: 0.0003, baseMean: 3421 },
  { gene: "BRCA1", lfc: -2.5, padj: 0.0002, baseMean: 4123 },
  { gene: "EGFR", lfc: 2.1, padj: 0.002, baseMean: 2876 },
  { gene: "KRAS", lfc: -1.9, padj: 0.004, baseMean: 3654 },
  { gene: "PIK3CA", lfc: 1.8, padj: 0.006, baseMean: 2234 },
  { gene: "PTEN", lfc: -2.2, padj: 0.001, baseMean: 3987 },
  { gene: "BRAF", lfc: 1.6, padj: 0.008, baseMean: 1987 },
  { gene: "AKT1", lfc: 2.4, padj: 0.0005, baseMean: 4321 },
  { gene: "NOTCH1", lfc: -1.7, padj: 0.007, baseMean: 2543 },
  { gene: "MAPK1", lfc: 1.5, padj: 0.01, baseMean: 3123 },
  { gene: "STAT3", lfc: 2.6, padj: 0.0004, baseMean: 3876 },
  { gene: "JAK2", lfc: -1.8, padj: 0.005, baseMean: 2765 },
  { gene: "VEGFA", lfc: 2.9, padj: 0.0002, baseMean: 4567 },
  { gene: "HIF1A", lfc: -2.1, padj: 0.003, baseMean: 3234 },
];

const pcaData = [
  { sample: "Control_1", pc1: -15, pc2: 8, condition: "control" },
  { sample: "Control_2", pc1: -18, pc2: 6, condition: "control" },
  { sample: "Control_3", pc1: -16, pc2: 10, condition: "control" },
  { sample: "Treated_1", pc1: 14, pc2: -7, condition: "treated" },
  { sample: "Treated_2", pc1: 16, pc2: -9, condition: "treated" },
  { sample: "Treated_3", pc1: 15, pc2: -6, condition: "treated" },
];

const stats = {
  totalGenes: 18432,
  upregulated: 1234,
  downregulated: 987,
  notSignificant: 16211,
};

export default function OmicsDemo() {
  const [padjThreshold, setPadjThreshold] = useState(0.05);
  const [lfcThreshold, setLfcThreshold] = useState(1.0);
  const [selectedTab, setSelectedTab] = useState<"volcano" | "pca" | "table">("volcano");

  const filteredData = sampleData.filter(
    (d) => d.padj < padjThreshold && Math.abs(d.lfc) > lfcThreshold
  );

  return (
    <div className="min-h-screen bg-gradient-to-br from-[#f0f4ff] via-white to-[#e6f0ff]">
      {/* Header */}
      <header className="border-b border-stone-200 bg-white/80 backdrop-blur">
        <div className="mx-auto max-w-7xl px-6 py-4">
          <div className="flex items-center justify-between">
            <Link
              href="/"
              className="text-sm font-semibold text-[var(--accent)] hover:text-[var(--accent-dark)]"
            >
              ← Back to ΕΚΦΑΝΣΙΣ
            </Link>
            <span className="text-xs uppercase tracking-[0.2em] text-stone-500">
              Interactive Demo
            </span>
          </div>
        </div>
      </header>

      {/* Hero Section */}
      <div className="mx-auto max-w-7xl px-6 py-12">
        <div className="text-center">
          <div className="inline-flex items-center gap-2 rounded-full bg-white px-4 py-2 text-xs font-semibold uppercase tracking-[0.2em] text-[var(--accent)] shadow-sm">
            <span>🧬</span>
            <span>Omics & Bioinformatics Demo</span>
          </div>
          <h1 className="mt-6 text-4xl font-semibold text-stone-900 sm:text-5xl">
            RNA-seq Differential Expression Analysis
          </h1>
          <p className="mx-auto mt-4 max-w-2xl text-lg text-stone-600">
            Interactive visualization of RNA-seq results. Explore differential expression data, adjust thresholds, and query results in real-time.
          </p>
        </div>

        {/* Stats Cards */}
        <div className="mt-12 grid gap-4 sm:grid-cols-2 lg:grid-cols-4">
          <div className="rounded-2xl border border-stone-200 bg-white p-6 shadow-sm">
            <p className="text-sm font-semibold uppercase tracking-[0.2em] text-stone-500">
              Total Genes
            </p>
            <p className="mt-2 text-3xl font-bold text-stone-900">
              {stats.totalGenes.toLocaleString()}
            </p>
          </div>
          <div className="rounded-2xl border border-stone-200 bg-white p-6 shadow-sm">
            <p className="text-sm font-semibold uppercase tracking-[0.2em] text-red-600">
              Upregulated
            </p>
            <p className="mt-2 text-3xl font-bold text-red-600">
              {stats.upregulated.toLocaleString()}
            </p>
          </div>
          <div className="rounded-2xl border border-stone-200 bg-white p-6 shadow-sm">
            <p className="text-sm font-semibold uppercase tracking-[0.2em] text-blue-600">
              Downregulated
            </p>
            <p className="mt-2 text-3xl font-bold text-blue-600">
              {stats.downregulated.toLocaleString()}
            </p>
          </div>
          <div className="rounded-2xl border border-stone-200 bg-white p-6 shadow-sm">
            <p className="text-sm font-semibold uppercase tracking-[0.2em] text-stone-500">
              Not Significant
            </p>
            <p className="mt-2 text-3xl font-bold text-stone-900">
              {stats.notSignificant.toLocaleString()}
            </p>
          </div>
        </div>

        {/* Controls */}
        <div className="mt-8 rounded-2xl border border-stone-200 bg-white p-6 shadow-sm">
          <p className="mb-4 text-sm font-semibold uppercase tracking-[0.2em] text-stone-700">
            Filter Controls
          </p>
          <div className="grid gap-6 sm:grid-cols-2">
            <div>
              <label className="mb-2 block text-sm font-medium text-stone-700">
                Adjusted P-value Threshold: {padjThreshold}
              </label>
              <input
                type="range"
                min="0.001"
                max="0.1"
                step="0.001"
                value={padjThreshold}
                onChange={(e) => setPadjThreshold(parseFloat(e.target.value))}
                className="w-full"
              />
              <div className="mt-1 flex justify-between text-xs text-stone-500">
                <span>0.001</span>
                <span>0.1</span>
              </div>
            </div>
            <div>
              <label className="mb-2 block text-sm font-medium text-stone-700">
                Log2 Fold Change Threshold: ±{lfcThreshold.toFixed(1)}
              </label>
              <input
                type="range"
                min="0.5"
                max="3.0"
                step="0.1"
                value={lfcThreshold}
                onChange={(e) => setLfcThreshold(parseFloat(e.target.value))}
                className="w-full"
              />
              <div className="mt-1 flex justify-between text-xs text-stone-500">
                <span>0.5</span>
                <span>3.0</span>
              </div>
            </div>
          </div>
          <div className="mt-4 rounded-lg bg-blue-50 p-3 text-sm text-blue-800">
            <strong>{filteredData.length}</strong> genes pass current thresholds
            (padj &lt; {padjThreshold.toFixed(3)}, |LFC| &gt; {lfcThreshold.toFixed(1)})
          </div>
        </div>

        {/* Tabs */}
        <div className="mt-8 flex gap-2 border-b border-stone-200">
          <button
            onClick={() => setSelectedTab("volcano")}
            className={`px-6 py-3 text-sm font-semibold transition ${
              selectedTab === "volcano"
                ? "border-b-2 border-[var(--accent)] text-[var(--accent)]"
                : "text-stone-500 hover:text-stone-700"
            }`}
          >
            Volcano Plot
          </button>
          <button
            onClick={() => setSelectedTab("pca")}
            className={`px-6 py-3 text-sm font-semibold transition ${
              selectedTab === "pca"
                ? "border-b-2 border-[var(--accent)] text-[var(--accent)]"
                : "text-stone-500 hover:text-stone-700"
            }`}
          >
            PCA Plot
          </button>
          <button
            onClick={() => setSelectedTab("table")}
            className={`px-6 py-3 text-sm font-semibold transition ${
              selectedTab === "table"
                ? "border-b-2 border-[var(--accent)] text-[var(--accent)]"
                : "text-stone-500 hover:text-stone-700"
            }`}
          >
            Data Table
          </button>
        </div>

        {/* Content */}
        <div className="mt-8">
          {selectedTab === "volcano" && (
            <VolcanoPlot data={sampleData} padjThreshold={padjThreshold} lfcThreshold={lfcThreshold} />
          )}
          {selectedTab === "pca" && <PCAPlot data={pcaData} />}
          {selectedTab === "table" && <DataTable data={filteredData} />}
        </div>

        {/* Technologies */}
        <div className="mt-12 rounded-2xl border border-stone-200 bg-white p-8 shadow-sm">
          <h2 className="text-2xl font-semibold text-stone-900">Technologies Used</h2>
          <div className="mt-6 grid gap-6 sm:grid-cols-3">
            <div>
              <h3 className="font-semibold text-[var(--accent)]">Workflow Management</h3>
              <ul className="mt-2 space-y-1 text-sm text-stone-600">
                <li>• Snakemake</li>
                <li>• Nextflow</li>
                <li>• Docker containers</li>
              </ul>
            </div>
            <div>
              <h3 className="font-semibold text-[var(--accent)]">Analysis Tools</h3>
              <ul className="mt-2 space-y-1 text-sm text-stone-600">
                <li>• Python + Pandas</li>
                <li>• SciPy + scikit-learn</li>
                <li>• DESeq2 (R)</li>
              </ul>
            </div>
            <div>
              <h3 className="font-semibold text-[var(--accent)]">Data Storage</h3>
              <ul className="mt-2 space-y-1 text-sm text-stone-600">
                <li>• PostgreSQL</li>
                <li>• DuckDB</li>
                <li>• RESTful APIs</li>
              </ul>
            </div>
          </div>
        </div>

        {/* CTA */}
        <div className="mt-12 rounded-2xl border border-[var(--accent)]/30 bg-gradient-to-br from-[#e0ebff] to-white p-8 text-center shadow-lg">
          <h2 className="text-2xl font-semibold text-stone-900">
            Interested in Custom Omics Analysis?
          </h2>
          <p className="mx-auto mt-4 max-w-2xl text-stone-600">
            We build reproducible bioinformatics pipelines, interactive dashboards, and data infrastructure for genomics, transcriptomics, proteomics, and more.
          </p>
          <div className="mt-6 flex justify-center gap-4">
            <a
              href="mailto:ekfansis@gmail.com?subject=Omics%20Services%20Inquiry"
              className="inline-flex items-center justify-center rounded-full bg-[var(--accent)] px-6 py-3 text-sm font-semibold text-white transition hover:bg-[var(--accent-dark)]"
            >
              Contact Us
            </a>
            <Link
              href="/"
              className="inline-flex items-center justify-center rounded-full border border-[var(--accent)] px-6 py-3 text-sm font-semibold text-[var(--accent)] transition hover:bg-[var(--accent)] hover:text-white"
            >
              Learn More
            </Link>
          </div>
        </div>
      </div>
    </div>
  );
}

function VolcanoPlot({
  data,
  padjThreshold,
  lfcThreshold,
}: {
  data: typeof sampleData;
  padjThreshold: number;
  lfcThreshold: number;
}) {
  const maxLogPadj = Math.max(...data.map((d) => -Math.log10(d.padj)));
  const maxLfc = Math.max(...data.map((d) => Math.abs(d.lfc)));

  return (
    <div className="rounded-2xl border border-stone-200 bg-white p-8 shadow-sm">
      <h3 className="mb-6 text-xl font-semibold text-stone-900">
        Volcano Plot: Treated vs Control
      </h3>
      <div className="relative" style={{ height: "500px" }}>
        <svg width="100%" height="100%" viewBox="0 0 800 500">
          {/* Axes */}
          <line x1="50" y1="450" x2="750" y2="450" stroke="#d4d4d4" strokeWidth="2" />
          <line x1="400" y1="50" x2="400" y2="450" stroke="#d4d4d4" strokeWidth="2" />

          {/* Threshold lines */}
          <line
            x1="50"
            y1={450 - (-Math.log10(padjThreshold) / maxLogPadj) * 400}
            x2="750"
            y2={450 - (-Math.log10(padjThreshold) / maxLogPadj) * 400}
            stroke="#9ca3af"
            strokeWidth="1"
            strokeDasharray="5,5"
          />
          <line
            x1={400 + (lfcThreshold / maxLfc) * 300}
            y1="50"
            x2={400 + (lfcThreshold / maxLfc) * 300}
            y2="450"
            stroke="#9ca3af"
            strokeWidth="1"
            strokeDasharray="5,5"
          />
          <line
            x1={400 - (lfcThreshold / maxLfc) * 300}
            y1="50"
            x2={400 - (lfcThreshold / maxLfc) * 300}
            y2="450"
            stroke="#9ca3af"
            strokeWidth="1"
            strokeDasharray="5,5"
          />

          {/* Data points */}
          {data.map((d, i) => {
            const x = 400 + (d.lfc / maxLfc) * 300;
            const y = 450 - (-Math.log10(d.padj) / maxLogPadj) * 400;
            const isSignificant = d.padj < padjThreshold && Math.abs(d.lfc) > lfcThreshold;
            const color = isSignificant
              ? d.lfc > 0
                ? "#ef4444"
                : "#3b82f6"
              : "#d4d4d4";

            return (
              <g key={i}>
                <circle cx={x} cy={y} r={isSignificant ? 6 : 4} fill={color} opacity={0.7} />
                {isSignificant && (
                  <text x={x + 10} y={y + 4} fontSize="10" fill="#4b5563">
                    {d.gene}
                  </text>
                )}
              </g>
            );
          })}

          {/* Axis labels */}
          <text x="400" y="480" textAnchor="middle" fontSize="14" fill="#6b7280">
            log2 Fold Change
          </text>
          <text
            x="20"
            y="250"
            textAnchor="middle"
            fontSize="14"
            fill="#6b7280"
            transform="rotate(-90, 20, 250)"
          >
            -log10(adjusted p-value)
          </text>
        </svg>
      </div>
    </div>
  );
}

function PCAPlot({ data }: { data: typeof pcaData }) {
  return (
    <div className="rounded-2xl border border-stone-200 bg-white p-8 shadow-sm">
      <h3 className="mb-6 text-xl font-semibold text-stone-900">
        PCA Plot: Sample Clustering
      </h3>
      <div className="relative" style={{ height: "500px" }}>
        <svg width="100%" height="100%" viewBox="0 0 800 500">
          {/* Axes */}
          <line x1="50" y1="250" x2="750" y2="250" stroke="#d4d4d4" strokeWidth="2" />
          <line x1="400" y1="50" x2="400" y2="450" stroke="#d4d4d4" strokeWidth="2" />

          {/* Data points */}
          {data.map((d, i) => {
            const x = 400 + d.pc1 * 15;
            const y = 250 - d.pc2 * 15;
            const color = d.condition === "control" ? "#3b82f6" : "#ef4444";

            return (
              <g key={i}>
                <circle cx={x} cy={y} r={10} fill={color} opacity={0.7} />
                <text x={x + 15} y={y + 4} fontSize="12" fill="#4b5563">
                  {d.sample}
                </text>
              </g>
            );
          })}

          {/* Legend */}
          <circle cx="650" cy="80" r="8" fill="#3b82f6" opacity={0.7} />
          <text x="665" y="85" fontSize="12" fill="#4b5563">
            Control
          </text>
          <circle cx="650" cy="110" r="8" fill="#ef4444" opacity={0.7} />
          <text x="665" y="115" fontSize="12" fill="#4b5563">
            Treated
          </text>

          {/* Axis labels */}
          <text x="400" y="480" textAnchor="middle" fontSize="14" fill="#6b7280">
            PC1 (42.3% variance)
          </text>
          <text
            x="20"
            y="250"
            textAnchor="middle"
            fontSize="14"
            fill="#6b7280"
            transform="rotate(-90, 20, 250)"
          >
            PC2 (28.7% variance)
          </text>
        </svg>
      </div>
    </div>
  );
}

function DataTable({ data }: { data: typeof sampleData }) {
  return (
    <div className="rounded-2xl border border-stone-200 bg-white shadow-sm">
      <div className="overflow-x-auto">
        <table className="w-full">
          <thead className="border-b border-stone-200 bg-stone-50">
            <tr>
              <th className="px-6 py-3 text-left text-xs font-semibold uppercase tracking-wider text-stone-600">
                Gene
              </th>
              <th className="px-6 py-3 text-left text-xs font-semibold uppercase tracking-wider text-stone-600">
                Log2 FC
              </th>
              <th className="px-6 py-3 text-left text-xs font-semibold uppercase tracking-wider text-stone-600">
                Adj P-value
              </th>
              <th className="px-6 py-3 text-left text-xs font-semibold uppercase tracking-wider text-stone-600">
                Base Mean
              </th>
              <th className="px-6 py-3 text-left text-xs font-semibold uppercase tracking-wider text-stone-600">
                Status
              </th>
            </tr>
          </thead>
          <tbody className="divide-y divide-stone-200">
            {data.map((row, i) => (
              <tr key={i} className="hover:bg-stone-50">
                <td className="whitespace-nowrap px-6 py-4 text-sm font-medium text-stone-900">
                  {row.gene}
                </td>
                <td
                  className={`whitespace-nowrap px-6 py-4 text-sm font-semibold ${
                    row.lfc > 0 ? "text-red-600" : "text-blue-600"
                  }`}
                >
                  {row.lfc > 0 ? "+" : ""}
                  {row.lfc.toFixed(2)}
                </td>
                <td className="whitespace-nowrap px-6 py-4 text-sm text-stone-600">
                  {row.padj.toExponential(2)}
                </td>
                <td className="whitespace-nowrap px-6 py-4 text-sm text-stone-600">
                  {row.baseMean.toLocaleString()}
                </td>
                <td className="whitespace-nowrap px-6 py-4">
                  <span
                    className={`inline-flex rounded-full px-2 py-1 text-xs font-semibold ${
                      row.lfc > 0
                        ? "bg-red-100 text-red-800"
                        : "bg-blue-100 text-blue-800"
                    }`}
                  >
                    {row.lfc > 0 ? "Upregulated" : "Downregulated"}
                  </span>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
