# Omics Demo - Live Interactive Visualization

## Online Demo

🎉 **Live URL:** `https://outwardly.net/omics-demo`

The demo is now accessible online! Πελάτες μπορούν να το δουν απευθείας στο browser τους.

## What's Included

### 1. **Interactive Visualizations**

#### Volcano Plot
- Real-time filtering με sliders
- Adjustable p-value threshold (0.001 - 0.1)
- Adjustable log2 fold change threshold (0.5 - 3.0)
- Color-coded genes (upregulated/downregulated)
- Gene labels για significant hits

#### PCA Plot
- Sample clustering visualization
- Control vs Treated groups
- PC1 και PC2 με explained variance
- Interactive legend

#### Data Table
- Sortable columns
- Gene names, log2FC, adjusted p-values
- Base mean expression values
- Status badges (upregulated/downregulated)

### 2. **Real-time Statistics**
- Total genes analyzed
- Number of upregulated genes
- Number of downregulated genes
- Non-significant genes

### 3. **Interactive Controls**
Δυναμικά sliders που επιτρέπουν στον χρήστη να:
- Αλλάξει το significance threshold
- Αλλάξει το fold change cutoff
- Δει πόσα genes περνούν τα criteria

### 4. **Sample Data**
15 cancer-related genes με realistic:
- Log2 fold changes
- Adjusted p-values
- Base mean expression
- TP53, MYC, BRCA1, EGFR, KRAS, κ.ά.

## Technical Features

### Built With
- **Next.js 15** - React framework
- **TypeScript** - Type safety
- **Tailwind CSS** - Styling
- **SVG** - Custom visualizations
- **Client-side rendering** - Fast, responsive

### Performance
- Static page generation
- No external dependencies για plots
- Lightweight (~6.64 KB route size)
- Mobile responsive

## How to Access

### From Main Page
1. Πήγαινε στο https://outwardly.net
2. Scroll στο "Services" section
3. Βρες την κάρτα "🧬 Omics & Bioinformatics"
4. Κάνε κλικ στο "View Interactive Demo →"

### Direct Link
https://outwardly.net/omics-demo

## Code Structure

```
src/app/omics-demo/
└── page.tsx              # Main demo page με όλα τα components
    ├── OmicsDemo         # Main component με state management
    ├── VolcanoPlot       # SVG volcano plot visualization
    ├── PCAPlot           # SVG PCA plot visualization
    └── DataTable         # Sortable data table
```

## Sample Data Format

```typescript
{
  gene: "TP53",           // Gene symbol
  lfc: 3.2,              // Log2 fold change
  padj: 0.0001,          // Adjusted p-value
  baseMean: 5432         // Base mean expression
}
```

## Future Enhancements

Potential additions για production demo:
- [ ] Real public dataset (GEO/TCGA)
- [ ] Download results as CSV
- [ ] Heatmap visualization
- [ ] Gene set enrichment preview
- [ ] API endpoint για queries
- [ ] User upload functionality
- [ ] Pathway analysis
- [ ] Integration με public databases (GeneCards, NCBI)

## Marketing Benefits

✅ **Tangible proof** - Δείχνει το product σε δράση
✅ **Interactive** - Πελάτες μπορούν να το δοκιμάσουν
✅ **Professional** - Polished UI/UX
✅ **Educational** - Εξηγεί τι κάνουμε
✅ **No friction** - Δεν χρειάζεται signup/login
✅ **Fast** - Loads instantly

## SEO & Discoverability

Το demo page έχει:
- Descriptive headings
- Keywords: genomics, RNA-seq, bioinformatics
- Clear CTAs
- Links στην κύρια σελίδα
- Contact information

## Analytics

Με το Google Analytics που προσθέσαμε, μπορείς να track:
- Page views στο demo
- Time on page
- Interaction με controls
- Conversion rate (demo → contact)

---

## Quick Start για Development

```bash
# Start dev server
npm run dev

# Visit demo
open http://localhost:3000/omics-demo

# Build for production
npm run build
npm start
```

## Deployment

Αυτό το demo είναι έτοιμο για production! Deploy με:
- **Vercel** (recommended για Next.js)
- **Netlify**
- **AWS Amplify**
- **Custom server**

---

**Created by OutWardly**
Contact: hello@outwardly.net
