"use client";

import Image from "next/image";
import { Suspense, useCallback, useEffect, useMemo, useState } from "react";
import { usePathname, useRouter, useSearchParams } from "next/navigation";
import { ProjectCarousel } from "./components/ProjectCarousel";
import { EUFundingBanner } from "./components/EUFundingBanner";

const EMAIL = "hello@ekfansis.com";
const PHONE = "+306945415350";

const translations = {
  el: {
    nav: {
      services: "Υπηρεσίες",
      projects: "Projects",
      methodology: "Μεθοδολογία",
      whyUs: "Γιατί εμάς",
      philosophy: "Φιλοσοφία",
      contact: "Επικοινωνία",
      menuOpen: "Άνοιγμα μενού",
      menuClose: "Κλείσιμο μενού",
    },
    hero: {
      heading: "Χτίζουμε custom ψηφιακές εμπειρίες και λογισμικό για την ομάδα σας.",
      body: "Από τη στρατηγική έως την υποστήριξη, η ομάδα μας σχεδιάζει, υλοποιεί και εξελίσσει την παρουσία σας στο διαδίκτυο. Μιλάμε τη γλώσσα της επιχείρησης και μεταφράζουμε τις ανάγκες σας σε κώδικα, design και μετρήσιμα αποτελέσματα.",
      primaryCta: "Κλείστε ένα ραντεβού",
      secondaryCta: "Δείτε τις υπηρεσίες μας",
    },
    summaryCard: {
      title: "Στόχος μας;",
      body: "Να γίνει το website σας ο καλύτερος πωλητής της εταιρίας. Επιλέγουμε τεχνολογίες που αντέχουν στον χρόνο, διασφαλίζουμε υψηλή απόδοση και επενδύουμε σε CX/UX που ξεχωρίζει.",
      bullets: [
        "Agile μεθοδολογία και εβδομαδιαία reports",
        "Integrations με τα εργαλεία που ήδη χρησιμοποιείτε",
        "Διαρκής παρακολούθηση και βελτιστοποιήσεις",
      ],
    },
    customSoftware: {
      label: "Custom Software Delivery",
      heading: "Από την ιδέα σε παραγωγή με end-to-end ανάπτυξη λογισμικού.",
      description:
        "Συνδυάζουμε product discovery, UX design και full-stack ανάπτυξη για να υλοποιήσουμε εργαλεία που λύνουν συγκεκριμένα επιχειρησιακά προβλήματα. Σχεδιάζουμε αρχιτεκτονική, στήνουμε pipelines και παραδίδουμε συνεχείς εκδόσεις που μετρούν πραγματικά KPIs.",
      bullets: [
        "MVPs, portals και integrations με ERP / CRM / third-party APIs.",
        "Τεχνικός σχεδιασμός, roadmaps ανά sprint και διαφάνεια σε κάθε release.",
        "DevOps + QA πρακτικές: CI/CD, automated testing και observability by default.",
      ],
      infoBlocks: [
        {
          title: "Tech Stack",
          body: "Next.js, Django (για μεγάλα eshop), Node.js, TypeScript, PostgreSQL, Prisma, Supabase, AWS, Vercel.",
        },
        {
          title: "Delivery",
          body: "Sprint-based υλοποίηση με εβδομαδιαία demos και shared dashboards για metrics.",
        },
        {
          title: "Συνέχεια",
          body: "SLA, υποστήριξη και R&D retainer ώστε το προϊόν να εξελίσσεται μαζί με την αγορά.",
        },
      ],
    },
    servicesSection: {
      label: "🧠 Ενότητα: Υπηρεσίες / Services",
      heading: "Ένα οικοσύστημα λύσεων για την ψηφιακή σας παρουσία.",
      phoneLabel: "+30 694 541 5350",
      cooperationCta: "Συνεργασία",
      mailtoSubject: "Νέο project για το portfolio",
      secondaryLink: "Θέλω να συζητήσουμε →",
    },
    services: [
      {
        emoji: "🚀",
        title: "Custom Λογισμικό",
        paragraphs: [
          "Σχεδιάζουμε και υλοποιούμε web εφαρμογές, portals και εσωτερικά συστήματα προσαρμοσμένα στις ροές της ομάδας σου.",
          "Ξεκινάμε με discovery workshops, service blueprints και τεχνική ανάλυση ώστε κάθε feature να απαντά σε πραγματική ανάγκη.",
          "Παραδίδουμε με agile iterations, CI/CD pipelines και τεχνική τεκμηρίωση που επιτρέπει στην ομάδα σου να εξελίσσει το προϊόν.",
        ],
        tagline: "→ MVPs, enterprise εργαλεία και αυτοματισμοί ραμμένοι στα μέτρα σου.",
      },
      {
        emoji: "🕸️",
        title: "Web Design & Development",
        paragraphs: [
          "Φτιάχνουμε καθαρές, γρήγορες, ουσιαστικές ιστοσελίδες.",
          "Χωρίς περιττά animations, χωρίς corporate φλυαρία.",
          "Χτίζουμε με Next.js και Tailwind, και για μεγάλα eshop αξιοποιούμε Django — για να φορτώνει γρήγορα, να δείχνει όμορφα και να λειτουργεί για τους ανθρώπους που τη χρησιμοποιούν.",
        ],
        tagline: "→ Websites που δεν χρειάζονται manual.",
      },
      {
        emoji: "🌐",
        title: "Web Hosting",
        paragraphs: [
          "Σύγχρονη υποδομή φιλοξενίας με real-time παρακολούθηση, αυτόματα back-ups και SLA που βασίζεται σε πραγματικές ανάγκες — όχι σε υποσχέσεις marketing.",
          "Μπορείς να κοιμάσαι ήσυχος· το site σου δεν το κάνει.",
        ],
        tagline: "→ Uptime, stability και ανθρώπινη υποστήριξη.",
      },
      {
        emoji: "🔧",
        title: "Maintenance & Support",
        paragraphs: [
          "Ενημερώσεις, έλεγχοι ασφαλείας, monitoring.",
          "Κρατάμε τα project ζωντανά — γιατί τίποτα δεν είναι “τελειωμένο” όταν είναι online.",
          "Δουλεύουμε αθόρυβα στο background, ώστε εσύ να εστιάζεις σε αυτό που έχει αξία.",
        ],
        tagline: "→ Τεχνική φροντίδα χωρίς corporate εισιτήρια υποστήριξης.",
      },
      {
        emoji: "🧩",
        title: "Custom Tools & Integrations",
        paragraphs: [
          "Μικρές web εφαρμογές, APIs, αυτοματισμοί και integrations με υπηρεσίες τρίτων.",
          "Απλοποιούμε τις διαδικασίες σου με λογική, όχι με buzzwords.",
        ],
        tagline: "→ Ό,τι χρειάζεσαι, χωρίς περιττά layers.",
      },
      {
        emoji: "🧬",
        title: "Omics & Bioinformatics",
        paragraphs: [
          "Ειδικευόμαστε σε genomics, transcriptomics, proteomics, metabolomics και epigenomics.",
          "Χτίζουμε reproducible workflows με Snakemake και Nextflow, αναλύουμε δημόσια datasets (RNA-seq, ChIP-seq, WGS) και εξάγουμε expression matrices.",
          "Στατιστική ανάλυση με Python/Pandas/SciPy/DESeq2, metadata σε PostgreSQL ή DuckDB, και custom dashboards για visualization.",
        ],
        tagline: "→ Data engineering + AI + web tech για τις omics επιστήμες.",
        demoLink: "/omics-demo",
      },
    ],
    projectsSection: {
      label: "Projects",
      heading: "Websites που ήδη ζουν εκεί έξω.",
      description:
        "Μερικές συνεργασίες που δείχνουν τι σημαίνει ΕΚΦΑΝΣΙΣ στην πράξη: από marketplace πλατφόρμες μέχρι media εμπειρίες με χιλιάδες χρήστες.",
      primaryButton: "Συνεργασία",
      secondaryLink: "Θέλω να συζητήσουμε →",
      mailtoSubject: "Νέο project για το portfolio",
      carouselPrev: "Προηγούμενο project",
      carouselNext: "Επόμενο project",
      carouselSlideLabel: "Project {current} από {total}",
    },
    projects: [
      {
        title: "dinalingerie.gr",
        subtitle: "E-commerce • Lingerie",
        description:
          "Χτίζουμε το νέο e-commerce για lingerie με καθαρή εμπειρία αγορών, omnichannel integrations και storytelling που αναδεικνύει το brand.",
        cta: "Σε εξέλιξη",
        href: "https://dinalingerie.gr",
        image: "/projects/dinalingerie.webp",
      },
      {
        title: "findteacher.gr",
        subtitle: "Marketplace • Εκπαίδευση",
        description:
          "Ξανασχεδιάσαμε τη διαδικασία εύρεσης καθηγητή: φίλτρα για επίπεδο, μάθημα και διαθεσιμότητα, με responsive UI που φορτώνει αστραπιαία και backend που χειρίζεται ασφαλείς κρατήσεις.",
        cta: "Δες το project →",
        href: "https://findteacher.gr",
        image: "/projects/findteacher.webp",
      },
      {
        title: "radioportal.me",
        subtitle: "Streaming • Media",
        description:
          "Ανανεώσαμε το radio streaming hub με custom player, real-time ενημέρωση προγράμματος και SEO-first αρχιτεκτονική ώστε να διαχειρίζεται χιλιάδες ακροατές χωρίς downtime.",
        cta: "Δες το project →",
        href: "https://radioportal.me",
        image: "/projects/radioportal.webp",
      },
      {
        title: "314project.gr",
        subtitle: "Events • Hospitality",
        description:
          "Δημιουργήσαμε την ψηφιακή παρουσία ενός πολυχώρου εκδηλώσεων, εκθέσεων και καφέ μπαρ με έμφαση στο storytelling, την εμπειρία του χρήστη και την εύκολη πλοήγηση.",
        cta: "Δες το project →",
        href: "https://314project.gr",
        image: "/projects/314project.webp",
      },
      {
        title: "kspkre.com",
        subtitle: "Real Estate • Marketplace",
        description:
          "Σχεδιάσαμε και αναπτύξαμε πλατφόρμα αγγελιών ακινήτων με advanced αναζήτηση, φίλτρα κατηγορίας και περιοχής, και καθαρό UI για αγοραστές και επαγγελματίες του κλάδου.",
        cta: "Δες το project →",
        href: "https://kspkre.com",
        image: "/projects/kspkre.webp",
      },
      {
        title: "chess.outwardly.net",
        subtitle: "PWA • Σκάκι",
        description:
          "Εργαλεία ανάλυσης και εκπαίδευσης σκακιού ως δωρεάν PWA: το προσθέτεις στο κινητό σου μέσω Chrome και έχεις PGN viewer και εκπαίδευση κινήσεων — μια δωρεάν λύση που έλειπε παγκόσμια.",
        cta: "Δες το project →",
        href: "https://chess.outwardly.net",
        image: "/projects/chess.webp",
      },
    ],
    methodology: {
      label: "Μεθοδολογία",
      heading: "Ένα συνεργατικό ταξίδι, με διαφάνεια και επίκεντρο τον χρήστη.",
      description:
        "Συνδυάζουμε στρατηγική, σχεδίαση και ανάπτυξη σε έναν κύκλο ζωής που προσαρμόζεται στις ανάγκες σας. Από την πρώτη συνάντηση έως το go-live και τα συνεχή iterations, είμαστε η ομάδα που θέλετε στο πλευρό σας.",
      steps: [
        {
          title: "Στρατηγική & Σχεδιασμός",
          text: "Καταγραφή στόχων, έρευνα κοινού και workshops για το brand και τα προϊόντα σας.",
        },
        {
          title: "Υλοποίηση & Ποιότητα",
          text: "Πλήρης ανάπτυξη με συνεχείς ελέγχους, αυτοματοποιημένα tests και διαφάνεια στην πρόοδο.",
          highlight: "QA • Monitoring • Green Deploys",
        },
        {
          title: "Λανσάρισμα & Βελτιστοποίηση",
          text: "Παράδοση, φιλοξενία και παρακολούθηση με βελτιώσεις βάσει δεδομένων και real-user metrics.",
        },
      ],
    },
    whyUs: {
      label: "Γιατί εμάς",
      heading: "Εμπιστευτείτε μια ομάδα που συνδυάζει δημιουργικότητα και τεχνογνωσία.",
      description:
        "Κάθε project συνοδεύεται από dedicated project manager, senior engineers και designers. Ερευνούμε, σχεδιάζουμε, υλοποιούμε και υποστηρίζουμε με πάθος για ποιότητα.",
      bullets: [
        {
          title: "Product & Engineering Squad",
          text: "Cross-functional ομάδα (product, design, dev, DevOps) που αναλαμβάνει discovery μέχρι rollout.",
        },
        {
          title: "Performance-first",
          text: "Lighthouse 90+, API benchmarks και observability dashboards με πραγματικά metrics.",
        },
        {
          title: "Μακροχρόνια σχέση",
          text: "Συμβόλαια υποστήριξης, SLA και roadmap sessions ανά τρίμηνο.",
        },
        {
          title: "Διαφανής κοστολόγηση",
          text: "Πακέτα και custom προσφορές με πλήρη ανάλυση ώρας και παραδοτέων.",
        },
      ],
    },
    philosophy: {
      label: "Φιλοσοφία",
      heading: "Μια μικρή ομάδα με βαθύ χρόνο για τα project της.",
      paragraphs: [
        "Η τεχνολογία δεν χρειάζεται να είναι απρόσωπη.",
        "Η ΕΚΦΑΝΣΙΣ είναι ένα μικρό, αυτόνομο creative studio.",
        "Χτίζουμε πράγματα που έχουν λόγο ύπαρξης — όχι απλώς “παρουσία στο web”.",
      ],
    },
    contact: {
      label: "Επικοινωνία",
      heading: "Ξεκινήστε σήμερα ένα project που θα μιλάει τη γλώσσα της αγοράς σας.",
      description:
        "Στείλτε μας μια σύντομη περιγραφή των στόχων σας και θα επιστρέψουμε με πρόταση, χρονοδιάγραμμα και ενδεικτικό προϋπολογισμό μέσα σε 2 εργάσιμες ημέρες.",
      introCta: "Κλείστε intro call →",
      form: {
        name: "Όνομα",
        email: "Email",
        company: "Εταιρεία (προαιρετικό)",
        message: "Περιγράψτε το project σας",
        submit: "Αποστολή μηνύματος",
        sending: "Αποστολή...",
        success: "Ευχαριστούμε! Θα επικοινωνήσουμε σύντομα.",
        error: "Κάτι πήγε στραβά. Δοκιμάστε ξανά ή στείλτε email.",
      },
    },
    footer: {
      taxId: "ΑΦΜ: 116201133",
      taxOffice: "ΔΟΥ: Πατρών",
      activityCodes: [
        "ΚΑΔ: 62.01.11.01 – Υπηρεσίες σχεδίασης και ανάπτυξης ιστοσελίδων",
        "ΚΑΔ: 63.11.11.00 – Υπηρεσίες φιλοξενίας ιστοσελίδων",
        "ΚΑΔ: 62.01.21.01 – Υπηρεσίες ανάπτυξης λογισμικού με βάση τις ανάγκες του πελάτη",
        "ΚΑΔ: 72.19.12.00 – Ερευνητικές και πειραματικές δραστηριότητες στη βιοτεχνολογία και βιοπληροφορική",
      ],
      contactTitle: "Επικοινωνία",
      phoneLabel: "+30 694 541 5350",
      introCta: "Κλείστε intro call →",
    },
  },
  en: {
    nav: {
      services: "Services",
      projects: "Projects",
      methodology: "Methodology",
      whyUs: "Why us",
      philosophy: "Philosophy",
      contact: "Contact",
      menuOpen: "Open menu",
      menuClose: "Close menu",
    },
    hero: {
      heading: "We craft custom digital experiences and software for your team.",
      body: "From strategy to support, our team designs, builds, and evolves your online presence. We speak the language of business and translate your needs into code, design, and measurable outcomes.",
      primaryCta: "Book a meeting",
      secondaryCta: "See our services",
    },
    summaryCard: {
      title: "Our focus",
      body: "To make your website your company’s best salesperson. We choose technologies that endure, ensure high performance, and invest in CX/UX that stands out.",
      bullets: [
        "Agile methodology and weekly reports",
        "Integrations with the tools you already use",
        "Continuous monitoring and optimisation",
      ],
    },
    customSoftware: {
      label: "Custom Software Delivery",
      heading: "From idea to production with end-to-end software development.",
      description:
        "We combine product discovery, UX design, and full-stack engineering to deliver tools that solve specific business problems. We design the architecture, set up the pipelines, and ship continuous releases that measure real KPIs.",
      bullets: [
        "MVPs, portals, and integrations with ERP / CRM / third-party APIs.",
        "Technical blueprinting, sprint roadmaps, and transparency in every release.",
        "DevOps + QA practices: CI/CD, automated testing, and observability by default.",
      ],
      infoBlocks: [
        {
          title: "Tech stack",
          body: "Next.js, Django (for large e-commerce), Node.js, TypeScript, PostgreSQL, Prisma, Supabase, AWS, Vercel.",
        },
        {
          title: "Delivery",
          body: "Sprint-based implementation with weekly demos and shared dashboards for metrics.",
        },
        {
          title: "Continuity",
          body: "SLA, support, and R&D retainers so the product evolves with the market.",
        },
      ],
    },
    servicesSection: {
      label: "🧠 Section: Services",
      heading: "An ecosystem of solutions for your digital presence.",
      phoneLabel: "+30 694 541 5350",
      cooperationCta: "Start a project",
      mailtoSubject: "New project for the portfolio",
      secondaryLink: "Let's talk →",
    },
    services: [
      {
        emoji: "🚀",
        title: "Custom Software",
        paragraphs: [
          "We design and deliver web applications, portals, and internal systems tailored to your team’s workflows.",
          "We begin with discovery workshops, service blueprints, and technical analysis so every feature answers a real need.",
          "We ship through agile iterations, CI/CD pipelines, and technical documentation that lets your team evolve the product.",
        ],
        tagline: "→ MVPs, enterprise tools, and automations built to measure.",
      },
      {
        emoji: "🕸️",
        title: "Web Design & Development",
        paragraphs: [
          "We build clean, fast, meaningful websites.",
          "No gratuitous animations, no corporate jargon.",
          "We use Next.js and Tailwind, and bring in Django for large e-commerce builds—so it loads quickly, looks great, and works for the people who use it.",
        ],
        tagline: "→ Websites that never need a manual.",
      },
      {
        emoji: "🌐",
        title: "Web Hosting",
        paragraphs: [
          "Modern hosting infrastructure with real-time monitoring, automatic backups, and an SLA based on actual needs—not marketing promises.",
          "Sleep easy; your site won’t.",
        ],
        tagline: "→ Uptime, stability, and human support.",
      },
      {
        emoji: "🔧",
        title: "Maintenance & Support",
        paragraphs: [
          "Updates, security checks, monitoring.",
          "We keep projects alive—because nothing online is ever “finished.”",
          "We work quietly in the background so you can focus on what matters.",
        ],
        tagline: "→ Technical care without corporate support tickets.",
      },
      {
        emoji: "🧩",
        title: "Custom Tools & Integrations",
        paragraphs: [
          "Lightweight web apps, APIs, automations, and third-party integrations.",
          "We simplify your processes with logic, not buzzwords.",
        ],
        tagline: "→ Exactly what you need, without extra layers.",
      },
      {
        emoji: "🧬",
        title: "Omics & Bioinformatics",
        paragraphs: [
          "We specialize in genomics, transcriptomics, proteomics, metabolomics, and epigenomics.",
          "We build reproducible workflows with Snakemake and Nextflow, analyze public datasets (RNA-seq, ChIP-seq, WGS), and extract expression matrices.",
          "Statistical analysis with Python/Pandas/SciPy/DESeq2, metadata in PostgreSQL or DuckDB, and custom dashboards for visualization.",
        ],
        tagline: "→ Data engineering + AI + web tech for omics sciences.",
        demoLink: "/omics-demo",
      },
    ],
    projectsSection: {
      label: "Projects",
      heading: "Websites already out in the wild.",
      description:
        "A few collaborations that show what EKFANSIS means in practice: from marketplace platforms to media experiences serving thousands of users.",
      primaryButton: "Start a project",
      secondaryLink: "Let's talk →",
      mailtoSubject: "New project for the portfolio",
      carouselPrev: "Previous project",
      carouselNext: "Next project",
      carouselSlideLabel: "Project {current} of {total}",
    },
    projects: [
      {
        title: "dinalingerie.gr",
        subtitle: "E-commerce • Lingerie",
        description:
          "We are crafting the new lingerie e-commerce experience with clean UX, omnichannel integrations, and storytelling that elevates the brand.",
        cta: "In progress",
        href: "https://dinalingerie.gr",
        image: "/projects/dinalingerie.webp",
      },
      {
        title: "findteacher.gr",
        subtitle: "Marketplace • Education",
        description:
          "We redesigned the teacher discovery journey: filters by level, subject, and availability, with a lightning-fast responsive UI and a backend that handles secure bookings.",
        cta: "View project →",
        href: "https://findteacher.gr",
        image: "/projects/findteacher.webp",
      },
      {
        title: "radioportal.me",
        subtitle: "Streaming • Media",
        description:
          "We refreshed the radio streaming hub with a custom player, real-time schedule updates, and SEO-first architecture that serves thousands of listeners without downtime.",
        cta: "View project →",
        href: "https://radioportal.me",
        image: "/projects/radioportal.webp",
      },
      {
        title: "314project.gr",
        subtitle: "Events • Hospitality",
        description:
          "We created the digital presence for a multifunctional event space, exhibition venue, and café bar with focus on storytelling, user experience, and seamless navigation.",
        cta: "View project →",
        href: "https://314project.gr",
        image: "/projects/314project.webp",
      },
      {
        title: "kspkre.com",
        subtitle: "Real Estate • Marketplace",
        description:
          "We designed and built a real estate listings platform with advanced search, category and location filters, and a clean UI serving both buyers and industry professionals.",
        cta: "View project →",
        href: "https://kspkre.com",
        image: "/projects/kspkre.webp",
      },
      {
        title: "chess.outwardly.net",
        subtitle: "PWA • Chess",
        description:
          "Free chess analysis and training tools as a PWA: add it to your phone via Chrome and get a PGN viewer and move training — a free solution that was missing worldwide.",
        cta: "View project →",
        href: "https://chess.outwardly.net",
        image: "/projects/chess.webp",
      },
    ],
    methodology: {
      label: "Methodology",
      heading: "A collaborative journey, transparent and user-focused.",
      description:
        "We combine strategy, design, and development in a lifecycle that adapts to your needs. From the first workshop to go-live and continuous iterations, we are the team you want by your side.",
      steps: [
        {
          title: "Strategy & Design",
          text: "Goal mapping, audience research, and workshops for your brand and products.",
        },
        {
          title: "Implementation & Quality",
          text: "End-to-end development with continuous reviews, automated tests, and transparent progress.",
          highlight: "QA • Monitoring • Green deploys",
        },
        {
          title: "Launch & Optimisation",
          text: "Delivery, hosting, and monitoring with improvements driven by data and real-user metrics.",
        },
      ],
    },
    whyUs: {
      label: "Why us",
      heading: "Trust a team that blends creativity with technical expertise.",
      description:
        "Every project includes a dedicated project manager, senior engineers, and designers. We research, design, build, and support with a passion for quality.",
      bullets: [
        {
          title: "Product & Engineering Squad",
          text: "A cross-functional team (product, design, dev, DevOps) that owns discovery through rollout.",
        },
        {
          title: "Performance-first",
          text: "Lighthouse 90+, API benchmarks, and observability dashboards anchored in real metrics.",
        },
        {
          title: "Long-term partnership",
          text: "Support retainers, SLAs, and quarterly roadmap sessions.",
        },
        {
          title: "Transparent costing",
          text: "Packages and custom proposals with a clear breakdown of hours and deliverables.",
        },
      ],
    },
    philosophy: {
      label: "Philosophy",
      heading: "A small team with deep focus for every project.",
      paragraphs: [
        "Technology doesn’t have to feel impersonal.",
        "EKFANSIS is a small, autonomous creative studio.",
        "We build things that have a reason to exist—not just a “web presence.”",
      ],
    },
    contact: {
      label: "Contact",
      heading: "Start a project that speaks your market's language.",
      description:
        "Send us a short outline of your goals and we'll reply with a proposal, timeline, and indicative budget within two business days.",
      introCta: "Schedule intro call →",
      form: {
        name: "Name",
        email: "Email",
        company: "Company (optional)",
        message: "Describe your project",
        submit: "Send message",
        sending: "Sending...",
        success: "Thank you! We'll be in touch soon.",
        error: "Something went wrong. Please try again or send an email.",
      },
    },
    footer: {
      taxId: "VAT ID: 116201133",
      taxOffice: "Tax Office: Patras",
      activityCodes: [
        "NACE 62.01.11.01 – Web design and development services",
        "NACE 63.11.11.00 – Web hosting services",
        "NACE 62.01.21.01 – Custom software development services",
        "NACE 72.19.12.00 – Research and experimental development in biotechnology and bioinformatics",
      ],
      contactTitle: "Contact",
      phoneLabel: "+30 694 541 5350",
      introCta: "Schedule intro call →",
    },
  },
} as const;

type Locale = keyof typeof translations;

const LOCALE_SWITCH: Record<Locale, { label: string; aria: string; next: Locale }> = {
  el: { label: "EN", aria: "Switch to English", next: "en" },
  en: { label: "ΕΛ", aria: "Switch to Greek", next: "el" },
};

export default function Home() {
  return (
    <Suspense fallback={<div className="min-h-screen bg-[var(--background)]" />}>
      <HomeContent />
    </Suspense>
  );
}

function HomeContent() {
  const searchParams = useSearchParams();
  const router = useRouter();
  const pathname = usePathname();

  const [locale, setLocaleState] = useState<Locale>(() =>
    searchParams.get("lang") === "en" ? "en" : "el",
  );
  const [menuOpen, setMenuOpen] = useState(false);

  useEffect(() => {
    const paramLocale = searchParams.get("lang") === "en" ? "en" : "el";
    setLocaleState(paramLocale);
  }, [searchParams]);

  const handleLocaleChange = useCallback(
    (nextLocale: Locale) => {
      setLocaleState(nextLocale);
      const params = new URLSearchParams(searchParams.toString());
      if (nextLocale === "el") {
        params.delete("lang");
      } else {
        params.set("lang", nextLocale);
      }
      const queryString = params.toString();
      router.replace(queryString ? `${pathname}?${queryString}` : pathname, {
        scroll: false,
      });
    },
    [pathname, router, searchParams],
  );

  const t = translations[locale];

  const introHref = useMemo(() => (locale === "en" ? "/intro?lang=en" : "/intro"), [locale]);

  return (
    <div className="min-h-screen bg-[var(--background)] text-[var(--foreground)]">
      {/* Sticky navigation */}
      <nav className="sticky top-0 z-50 border-b border-stone-900/10 bg-white/90 backdrop-blur" aria-label="Main navigation">
        <div className="mx-auto flex max-w-6xl items-center justify-between px-6 py-3 lg:px-8">
          <a href="#" className="text-sm font-bold tracking-wide text-stone-900">ΕΚΦΑΝΣΙΣ</a>

          {/* Desktop links */}
          <div className="hidden items-center gap-6 md:flex">
            {([
              ["#services", t.nav.services],
              ["#projects", t.nav.projects],
              ["#methodology", t.nav.methodology],
              ["#why-us", t.nav.whyUs],
              ["#philosophy", t.nav.philosophy],
              ["#contact", t.nav.contact],
            ] as const).map(([href, label]) => (
              <a
                key={href}
                href={href}
                className="text-sm text-stone-600 transition hover:text-[var(--accent)]"
              >
                {label}
              </a>
            ))}
            <button
              type="button"
              onClick={() => handleLocaleChange(LOCALE_SWITCH[locale].next)}
              className="rounded-full border border-[var(--accent)]/40 px-3 py-1 text-xs font-semibold uppercase tracking-[0.2em] text-[var(--accent)] transition hover:bg-[var(--accent)] hover:text-white"
              aria-label={LOCALE_SWITCH[locale].aria}
            >
              {LOCALE_SWITCH[locale].label}
            </button>
          </div>

          {/* Mobile: lang + hamburger */}
          <div className="flex items-center gap-3 md:hidden">
            <button
              type="button"
              onClick={() => handleLocaleChange(LOCALE_SWITCH[locale].next)}
              className="rounded-full border border-[var(--accent)]/40 px-3 py-1 text-xs font-semibold uppercase tracking-[0.2em] text-[var(--accent)]"
              aria-label={LOCALE_SWITCH[locale].aria}
            >
              {LOCALE_SWITCH[locale].label}
            </button>
            <button
              type="button"
              onClick={() => setMenuOpen(!menuOpen)}
              className="flex h-9 w-9 items-center justify-center rounded-lg text-stone-700 transition hover:bg-stone-100"
              aria-label={menuOpen ? t.nav.menuClose : t.nav.menuOpen}
              aria-expanded={menuOpen}
            >
              {menuOpen ? (
                <svg className="h-5 w-5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                  <path strokeLinecap="round" strokeLinejoin="round" d="M6 18L18 6M6 6l12 12" />
                </svg>
              ) : (
                <svg className="h-5 w-5" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth={2}>
                  <path strokeLinecap="round" strokeLinejoin="round" d="M4 6h16M4 12h16M4 18h16" />
                </svg>
              )}
            </button>
          </div>
        </div>

        {/* Mobile dropdown */}
        {menuOpen && (
          <div className="border-t border-stone-900/10 bg-white px-6 pb-4 pt-2 md:hidden">
            {([
              ["#services", t.nav.services],
              ["#projects", t.nav.projects],
              ["#methodology", t.nav.methodology],
              ["#why-us", t.nav.whyUs],
              ["#philosophy", t.nav.philosophy],
              ["#contact", t.nav.contact],
            ] as const).map(([href, label]) => (
              <a
                key={href}
                href={href}
                onClick={() => setMenuOpen(false)}
                className="block py-2 text-sm text-stone-600 transition hover:text-[var(--accent)]"
              >
                {label}
              </a>
            ))}
          </div>
        )}
      </nav>

      <header className="relative isolate overflow-hidden bg-white/80">
        <div className="absolute inset-0 -z-10 bg-[radial-gradient(circle_at_top,_rgba(37,99,235,0.22)_0%,_transparent_58%)]" />
        <div className="px-4 pt-3">
          <EUFundingBanner />
        </div>
        <div className="mx-auto max-w-6xl px-6 pb-24 pt-16 sm:pb-28 sm:pt-20 lg:px-8">
          <div className="mb-12">
            <Image
              src="/ekfansis-logo.svg"
              alt="ΕΚΦΑΝΣΙΣ logo"
              width={800}
              height={400}
              priority
              className="mx-auto w-full max-w-5xl"
            />
          </div>
          <div className="flex flex-col gap-12 lg:flex-row lg:items-center lg:gap-20">
          <div className="flex-1 space-y-6">
            <h1 className="text-4xl font-semibold leading-tight sm:text-5xl">{t.hero.heading}</h1>
            <p className="max-w-2xl text-base text-stone-700 sm:text-lg">{t.hero.body}</p>
            <div className="flex flex-col gap-3 sm:flex-row">
              <a
                className="inline-flex items-center justify-center rounded-full bg-[var(--accent)] px-6 py-3 text-sm font-semibold text-white transition hover:bg-[var(--accent-dark)]"
                href={`mailto:${EMAIL}`}
              >
                {t.hero.primaryCta}
              </a>
              <a
                className="inline-flex items-center justify-center rounded-full border border-[var(--accent)] px-6 py-3 text-sm font-semibold text-[var(--accent)] transition hover:bg-[var(--accent)] hover:text-white"
                href="#services"
              >
                {t.hero.secondaryCta}
              </a>
            </div>
          </div>
          <div className="flex max-w-md flex-1 flex-col gap-4 rounded-3xl border border-stone-900/10 bg-white p-6 shadow-lg backdrop-blur">
            <h2 className="text-lg font-semibold text-stone-700">{t.summaryCard.title}</h2>
            <p className="text-sm text-stone-700">{t.summaryCard.body}</p>
            <ul className="space-y-2 text-sm text-stone-600">
              {t.summaryCard.bullets.map((bullet) => (
                <li key={bullet}>• {bullet}</li>
              ))}
            </ul>
          </div>
          </div>
        </div>
      </header>

      <main id="main-content" role="main" className="mx-auto flex max-w-6xl flex-col gap-24 px-6 py-16 sm:px-8">
        <section className="grid gap-8 rounded-3xl border border-[var(--accent)]/30 bg-gradient-to-br from-white via-[#eef3ff] to-[#dbe7ff] p-10 shadow-lg lg:grid-cols-[1.15fr_0.85fr] lg:items-center">
          <div className="space-y-6">
            <p className="inline-flex items-center gap-2 text-xs font-semibold uppercase tracking-[0.3em] text-[var(--accent)]">
              {t.customSoftware.label}
            </p>
            <h2 className="text-3xl font-semibold sm:text-4xl">{t.customSoftware.heading}</h2>
            <p className="text-base text-stone-600 sm:text-lg">{t.customSoftware.description}</p>
            <ul className="space-y-3 text-sm text-stone-700">
              {t.customSoftware.bullets.map((item) => (
                <li key={item} className="flex items-start gap-3">
                  <span className="mt-1 h-2 w-2 rounded-full bg-[var(--accent)]" aria-hidden="true" />
                  {item}
                </li>
              ))}
            </ul>
          </div>
          <div className="grid gap-4 rounded-3xl border border-white/60 bg-white/80 p-6 shadow-lg backdrop-blur">
            {t.customSoftware.infoBlocks.map((block) => (
              <div key={block.title} className="space-y-2 border-t border-stone-200/70 pt-4 first:border-none first:pt-0">
                <p className="text-xs font-semibold uppercase tracking-[0.28em] text-[var(--accent)]">{block.title}</p>
                <p className="text-sm text-stone-600">{block.body}</p>
              </div>
            ))}
          </div>
        </section>

        <section id="services" className="space-y-10">
          <div className="flex flex-col gap-4 sm:flex-row sm:items-end sm:justify-between">
            <div>
              <p className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">{t.servicesSection.label}</p>
              <h2 className="text-3xl font-semibold sm:text-4xl">{t.servicesSection.heading}</h2>
            </div>
            <a
              className="text-sm font-semibold text-[var(--accent)] hover:text-[var(--accent-dark)]"
              href={`tel:${PHONE}`}
            >
              {t.servicesSection.phoneLabel}
            </a>
          </div>
          <div className="grid gap-6 sm:grid-cols-2">
            {t.services.map((service) => (
              <div
                key={service.title}
                className="group flex h-full flex-col gap-5 rounded-3xl border border-stone-900/10 bg-white p-6 shadow-sm transition hover:-translate-y-1 hover:border-[var(--accent)]/50 hover:shadow-lg"
              >
                <div>
                  <span className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">
                    {service.emoji}
                  </span>
                  <h3 className="mt-2 text-2xl font-semibold text-stone-900">{service.title}</h3>
                </div>
                <div className="space-y-3 text-sm text-stone-600">
                  {service.paragraphs.map((paragraph) => (
                    <p key={paragraph}>{paragraph}</p>
                  ))}
                </div>
                <p className="text-sm font-semibold text-[var(--accent)]">{service.tagline}</p>
                {"demoLink" in service && service.demoLink && (
                  <a
                    href={service.demoLink}
                    className="mt-auto inline-flex items-center gap-2 text-sm font-semibold text-[var(--accent)] transition hover:text-[var(--accent-dark)]"
                  >
                    View Interactive Demo →
                  </a>
                )}
              </div>
            ))}
          </div>
        </section>

        <section id="projects" className="space-y-10 rounded-3xl border border-stone-900/10 bg-white p-10 shadow-lg">
          <div className="flex flex-col gap-4 sm:flex-row sm:items-end sm:justify-between">
            <div className="space-y-2">
              <p className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">{t.projectsSection.label}</p>
              <h2 className="text-3xl font-semibold sm:text-4xl">{t.projectsSection.heading}</h2>
              <p className="text-base text-stone-600">{t.projectsSection.description}</p>
            </div>
            <div className="flex gap-3">
              <a
                className="inline-flex items-center justify-center rounded-full border border-[var(--accent)] px-5 py-2 text-xs font-semibold uppercase tracking-[0.2em] text-[var(--accent)] transition hover:bg-[var(--accent)] hover:text-white"
                href="#contact"
              >
                {t.projectsSection.primaryButton}
              </a>
              <a
                className="text-sm font-semibold text-[var(--accent)] hover:text-[var(--accent-dark)]"
                href={`mailto:${EMAIL}?subject=${encodeURIComponent(t.projectsSection.mailtoSubject)}`}
              >
                {t.projectsSection.secondaryLink}
              </a>
            </div>
          </div>
          <ProjectCarousel
            projects={t.projects}
            ariaLabels={{
              prev: t.projectsSection.carouselPrev,
              next: t.projectsSection.carouselNext,
              slideLabel: t.projectsSection.carouselSlideLabel,
            }}
          />
        </section>

        <section id="methodology" className="grid gap-10 lg:grid-cols-[1.1fr_0.9fr] lg:items-center">
          <div className="space-y-6">
            <p className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">{t.methodology.label}</p>
            <h2 className="text-3xl font-semibold sm:text-4xl">{t.methodology.heading}</h2>
            <p className="text-base text-stone-600">{t.methodology.description}</p>
          </div>
          <div className="space-y-6 rounded-3xl border border-stone-900/10 bg-white p-8 shadow-sm">
            {t.methodology.steps.map((step, index) => (
              <div key={step.title} className="space-y-2 border-l-2 border-[var(--accent)]/40 pl-6">
                <span className="text-xs font-semibold uppercase tracking-[0.3em] text-[var(--accent)]">
                  {locale === "el" ? `Βήμα ${index + 1}` : `Step ${index + 1}`}
                </span>
                <h3 className="text-lg font-semibold text-stone-800">{step.title}</h3>
                <p className="text-sm text-stone-600">{step.text}</p>
                {"highlight" in step && step.highlight ? (
                  <p className="text-xs font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">{step.highlight}</p>
                ) : null}
              </div>
            ))}
          </div>
        </section>

        <section id="why-us" className="grid gap-10 rounded-3xl border border-stone-900/10 bg-white p-10 shadow-sm lg:grid-cols-[1.1fr_0.9fr] lg:items-center">
          <div className="space-y-6">
            <p className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">{t.whyUs.label}</p>
            <h2 className="text-3xl font-semibold sm:text-4xl">{t.whyUs.heading}</h2>
            <p className="text-base text-stone-600">{t.whyUs.description}</p>
          </div>
          <ul className="space-y-4 text-sm text-stone-700">
            {t.whyUs.bullets.map((item) => (
              <li key={item.title} className="flex items-start gap-3 rounded-2xl border border-stone-900/10 bg-[#eef3ff] p-4">
                <span className="mt-0.5 text-lg text-[var(--accent)]">◆</span>
                <div>
                  <p className="font-semibold text-stone-900">{item.title}</p>
                  <p className="text-stone-600">{item.text}</p>
                </div>
              </li>
            ))}
          </ul>
        </section>

        <section id="philosophy" className="rounded-3xl border border-stone-900/10 bg-white p-10 shadow-sm">
          <div className="space-y-4">
            <p className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">
              {t.philosophy.label}
            </p>
            <h2 className="text-3xl font-semibold text-stone-900 sm:text-4xl">{t.philosophy.heading}</h2>
            <div className="space-y-3 text-base text-stone-600">
              {t.philosophy.paragraphs.map((paragraph) => (
                <p key={paragraph}>{paragraph}</p>
              ))}
            </div>
          </div>
        </section>

        <section
          id="contact"
          aria-labelledby="contact-heading"
          className="relative overflow-hidden rounded-3xl border border-[var(--accent)]/30 bg-gradient-to-br from-[#e0ebff] via-[#f3f7ff] to-[var(--background)] p-10"
        >
          <div className="absolute inset-0 -z-10 bg-[radial-gradient(circle_at_bottom_right,_rgba(37,99,235,0.18)_0%,_transparent_55%)]" aria-hidden="true" />
          <div className="space-y-6">
            <p className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">{t.contact.label}</p>
            <h2 id="contact-heading" className="text-3xl font-semibold sm:text-4xl">{t.contact.heading}</h2>
            <p className="max-w-2xl text-base text-stone-700 sm:text-lg">{t.contact.description}</p>
            <div className="flex flex-col gap-3 sm:flex-row">
              <a
                className="inline-flex items-center justify-center rounded-full bg-[var(--accent)] px-6 py-3 text-sm font-semibold text-white transition hover:bg-[var(--accent-dark)]"
                href={`mailto:${EMAIL}`}
                aria-label={`Αποστολή email στο ${EMAIL}`}
              >
                {EMAIL}
              </a>
            </div>
          </div>
        </section>
      </main>

      <footer role="contentinfo" aria-label="Στοιχεία επικοινωνίας" className="border-t border-stone-300 bg-[#e7edff]">
        <div className="mx-auto grid max-w-6xl gap-6 px-6 py-8 text-sm text-stone-600 sm:grid-cols-3 sm:px-8">
          <div className="space-y-1 text-stone-700">
            <p className="font-semibold text-stone-900">© 2025 ΕΚΦΑΝΣΙΣ</p>
            <p>{t.footer.taxId}</p>
            <p>{t.footer.taxOffice}</p>
          </div>
          <div className="space-y-1" aria-label="Κωδικοί δραστηριότητας">
            {t.footer.activityCodes.map((code) => (
              <p key={code}>{code}</p>
            ))}
          </div>
          <div className="space-y-1">
            <p className="font-semibold text-stone-900">{t.footer.contactTitle}</p>
            <a
              className="block hover:text-[var(--accent)]"
              href={`tel:${PHONE}`}
              aria-label={`Τηλεφωνική επικοινωνία: ${t.footer.phoneLabel}`}
            >
              {t.footer.phoneLabel}
            </a>
            <a
              className="block hover:text-[var(--accent)]"
              href={`mailto:${EMAIL}`}
              aria-label={`Αποστολή email στο ${EMAIL}`}
            >
              {EMAIL}
            </a>
          </div>
        </div>
      </footer>
    </div>
  );
}

interface ContactFormProps {
  translations: {
    name: string;
    email: string;
    company: string;
    message: string;
    submit: string;
    sending: string;
    success: string;
    error: string;
  };
}

function ContactForm({ translations: t }: ContactFormProps) {
  const [status, setStatus] = useState<"idle" | "sending" | "success" | "error">("idle");

  const handleSubmit = async (e: React.FormEvent<HTMLFormElement>) => {
    e.preventDefault();
    setStatus("sending");

    const form = e.currentTarget;
    const formData = new FormData(form);

    try {
      const response = await fetch("https://formspree.io/f/xzzeddpn", {
        method: "POST",
        body: formData,
        headers: {
          Accept: "application/json",
        },
      });

      if (response.ok) {
        setStatus("success");
        form.reset();
      } else {
        setStatus("error");
      }
    } catch {
      setStatus("error");
    }
  };

  return (
    <form
      onSubmit={handleSubmit}
      className="space-y-5 rounded-2xl border border-white/60 bg-white/90 p-6 shadow-lg backdrop-blur"
    >
      <div className="space-y-4">
        <div>
          <label htmlFor="name" className="mb-1.5 block text-sm font-medium text-stone-700">
            {t.name}
          </label>
          <input
            type="text"
            id="name"
            name="name"
            required
            className="w-full rounded-xl border border-stone-200 bg-white px-4 py-3 text-sm text-stone-900 transition placeholder:text-stone-400 focus:border-[var(--accent)] focus:outline-none focus:ring-2 focus:ring-[var(--accent)]/20"
          />
        </div>
        <div>
          <label htmlFor="email" className="mb-1.5 block text-sm font-medium text-stone-700">
            {t.email}
          </label>
          <input
            type="email"
            id="email"
            name="email"
            required
            className="w-full rounded-xl border border-stone-200 bg-white px-4 py-3 text-sm text-stone-900 transition placeholder:text-stone-400 focus:border-[var(--accent)] focus:outline-none focus:ring-2 focus:ring-[var(--accent)]/20"
          />
        </div>
        <div>
          <label htmlFor="company" className="mb-1.5 block text-sm font-medium text-stone-700">
            {t.company}
          </label>
          <input
            type="text"
            id="company"
            name="company"
            className="w-full rounded-xl border border-stone-200 bg-white px-4 py-3 text-sm text-stone-900 transition placeholder:text-stone-400 focus:border-[var(--accent)] focus:outline-none focus:ring-2 focus:ring-[var(--accent)]/20"
          />
        </div>
        <div>
          <label htmlFor="message" className="mb-1.5 block text-sm font-medium text-stone-700">
            {t.message}
          </label>
          <textarea
            id="message"
            name="message"
            rows={4}
            required
            className="w-full resize-none rounded-xl border border-stone-200 bg-white px-4 py-3 text-sm text-stone-900 transition placeholder:text-stone-400 focus:border-[var(--accent)] focus:outline-none focus:ring-2 focus:ring-[var(--accent)]/20"
          />
        </div>
      </div>

      {status === "success" && (
        <p className="rounded-xl bg-emerald-50 px-4 py-3 text-sm font-medium text-emerald-700">
          {t.success}
        </p>
      )}
      {status === "error" && (
        <p className="rounded-xl bg-red-50 px-4 py-3 text-sm font-medium text-red-700">
          {t.error}
        </p>
      )}

      <button
        type="submit"
        disabled={status === "sending"}
        className="w-full rounded-full bg-[var(--accent)] px-6 py-3 text-sm font-semibold text-white transition hover:bg-[var(--accent-dark)] disabled:cursor-not-allowed disabled:opacity-60"
      >
        {status === "sending" ? t.sending : t.submit}
      </button>
    </form>
  );
}
