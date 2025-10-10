import Image from "next/image";

export default function Home() {
  const services = [
    {
      emoji: "🕸️",
      title: "Web Design & Development",
      paragraphs: [
        "Φτιάχνουμε καθαρές, γρήγορες, ουσιαστικές ιστοσελίδες.",
        "Χωρίς περιττά animations, χωρίς corporate φλυαρία.",
        "Χτίζουμε με Next.js, Tailwind, και προσοχή στη λεπτομέρεια — για να φορτώνει γρήγορα, να δείχνει όμορφα και να λειτουργεί για τους ανθρώπους που τη χρησιμοποιούν.",
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
  ];

  const philosophy = {
    emoji: "🧭",
    title: "Φιλοσοφία",
    paragraphs: [
      "Η τεχνολογία δεν χρειάζεται να είναι απρόσωπη.",
      "Η Outwardly είναι ένας μικρός, αυτόνομος κόμβος δημιουργίας.",
      "Χτίζουμε πράγματα που έχουν λόγο ύπαρξης — όχι απλώς “παρουσία στο web”.",
    ],
  };

  const projects = [
    {
      title: "findteacher.gr",
      subtitle: "Marketplace • Εκπαίδευση",
      description:
        "Ξανασχεδιάσαμε τη διαδικασία εύρεσης καθηγητή: φίλτρα για επίπεδο, μάθημα και διαθεσιμότητα, με responsive UI που φορτώνει αστραπιαία και backend που χειρίζεται ασφαλείς κρατήσεις.",
      href: "https://findteacher.gr",
    },
    {
      title: "radioportal.me",
      subtitle: "Streaming • Media",
      description:
        "Ανανεώσαμε το radio streaming hub με custom player, real-time ενημέρωση προγράμματος και SEO-first αρχιτεκτονική ώστε να διαχειρίζεται χιλιάδες ακροατές χωρίς downtime.",
      href: "https://radioportal.me",
    },
  ];

  const milestones = [
    {
      title: "Στρατηγική & Σχεδιασμός",
      text: "Καταγραφή στόχων, έρευνα κοινού και workshops για το brand και τα προϊόντα σας.",
    },
    {
      title: "Υλοποίηση & Ποιότητα",
      text: "Πλήρης ανάπτυξη με συνεχείς ελέγχους, αυτοματοποιημένα tests και διαφάνεια στην πρόοδο.",
    },
    {
      title: "Λανσάρισμα & Βελτιστοποίηση",
      text: "Παράδοση, φιλοξενία και παρακολούθηση με βελτιώσεις βάσει δεδομένων και real-user metrics.",
    },
  ];

  return (
    <div className="min-h-screen bg-[var(--background)] text-[var(--foreground)]">
      <header className="relative isolate overflow-hidden bg-white/80">
        <div className="absolute inset-0 -z-10 bg-[radial-gradient(circle_at_top,_rgba(37,99,235,0.22)_0%,_transparent_58%)]" />
        <div className="mx-auto flex max-w-6xl flex-col gap-12 px-6 pb-24 pt-16 sm:pb-28 sm:pt-24 lg:flex-row lg:items-center lg:gap-20 lg:px-8">
          <div className="flex-1 space-y-6">
            <div>
              <Image
                src="/OutWardly_logo.png"
                alt="OutWardly logo"
                width={1024}
                height={1024}
                priority
                className="max-w-full rounded-xl shadow-lg"
              />
            </div>
            <h1 className="text-4xl font-semibold leading-tight sm:text-5xl">
              Χτίζουμε αξιόπιστες ψηφιακές εμπειρίες για την πρώτη σας εταιρία.
            </h1>
            <p className="max-w-2xl text-base text-stone-700 sm:text-lg">
              Από τη στρατηγική έως την υποστήριξη, η ομάδα μας σχεδιάζει,
              υλοποιεί και εξελίσσει την παρουσία σας στο διαδίκτυο. Μιλάμε τη
              γλώσσα της επιχείρησης και μεταφράζουμε τις ανάγκες σας σε
              κώδικα, design και μετρήσιμα αποτελέσματα.
            </p>
            <div className="flex flex-col gap-3 sm:flex-row">
              <a
                className="inline-flex items-center justify-center rounded-full bg-[var(--accent)] px-6 py-3 text-sm font-semibold text-white transition hover:bg-[var(--accent-dark)]"
                href="mailto:hello@outwardly.gr"
              >
                Κλείστε ένα ραντεβού
              </a>
              <a
                className="inline-flex items-center justify-center rounded-full border border-[var(--accent)] px-6 py-3 text-sm font-semibold text-[var(--accent)] transition hover:bg-[var(--accent)] hover:text-white"
                href="#services"
              >
                Δείτε τις υπηρεσίες μας
              </a>
            </div>
          </div>
          <div className="flex max-w-md flex-1 flex-col gap-4 rounded-3xl border border-stone-900/10 bg-white p-6 shadow-lg backdrop-blur">
            <h2 className="text-lg font-semibold text-stone-700">
              Στόχος μας;
            </h2>
            <p className="text-sm text-stone-700">
              Να γίνει το website σας ο καλύτερος πωλητής της εταιρίας. Επιλέγουμε
              τεχνολογίες που αντέχουν στον χρόνο, διασφαλίζουμε υψηλή απόδοση και
              επενδύουμε σε CX/UX που ξεχωρίζει.
            </p>
            <ul className="space-y-2 text-sm text-stone-600">
              <li>• Agile μεθοδολογία και εβδομαδιαία reports</li>
              <li>• Integrations με τα εργαλεία που ήδη χρησιμοποιείτε</li>
              <li>• Διαρκής παρακολούθηση και βελτιστοποιήσεις</li>
            </ul>
          </div>
        </div>
      </header>

      <main className="mx-auto flex max-w-6xl flex-col gap-24 px-6 py-16 sm:px-8">
        <section id="services" className="space-y-10">
          <div className="flex flex-col gap-4 sm:flex-row sm:items-end sm:justify-between">
            <div>
              <p className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">
                🧠 Ενότητα: Υπηρεσίες / Services
              </p>
              <h2 className="text-3xl font-semibold sm:text-4xl">
                Ένα οικοσύστημα λύσεων για την ψηφιακή σας παρουσία.
              </h2>
            </div>
            <a
              className="text-sm font-semibold text-[var(--accent)] hover:text-[var(--accent-dark)]"
              href="tel:+306945415350"
            >
              +30 694 541 5350
            </a>
          </div>
          <div className="grid gap-6 sm:grid-cols-2">
            {services.map((service) => (
              <div
                key={service.title}
                className="group flex h-full flex-col gap-5 rounded-3xl border border-stone-900/10 bg-white p-6 shadow-sm transition hover:-translate-y-1 hover:border-[var(--accent)]/50 hover:shadow-lg"
              >
                <div>
                  <span className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">
                    {service.emoji}
                  </span>
                  <h3 className="mt-2 text-2xl font-semibold text-stone-900">
                    {service.title}
                  </h3>
                </div>
                <div className="space-y-3 text-sm text-stone-600">
                  {service.paragraphs.map((paragraph) => (
                    <p key={paragraph}>{paragraph}</p>
                  ))}
                </div>
                <p className="text-sm font-semibold text-[var(--accent)]">{service.tagline}</p>
              </div>
            ))}
          </div>
        </section>

        <section className="space-y-10 rounded-3xl border border-stone-900/10 bg-white p-10 shadow-lg">
          <div className="flex flex-col gap-4 sm:flex-row sm:items-end sm:justify-between">
            <div className="space-y-2">
              <p className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">
                Projects
              </p>
              <h2 className="text-3xl font-semibold sm:text-4xl">
                Websites που ήδη ζουν εκεί έξω.
              </h2>
              <p className="text-base text-stone-600">
                Δύο συνεργασίες που δείχνουν τι σημαίνει Outwardly στην πράξη: από marketplace
                πλατφόρμες μέχρι media εμπειρίες με χιλιάδες χρήστες.
              </p>
            </div>
            <div className="flex gap-3">
              <a
                className="inline-flex items-center justify-center rounded-full border border-[var(--accent)] px-5 py-2 text-xs font-semibold uppercase tracking-[0.2em] text-[var(--accent)] transition hover:bg-[var(--accent)] hover:text-white"
                href="#contact"
              >
                Συνεργασία
              </a>
              <a
                className="text-sm font-semibold text-[var(--accent)] hover:text-[var(--accent-dark)]"
                href="mailto:hello@outwardly.gr?subject=Νέο%20project%20για%20το%20portfolio"
              >
                Θέλω να συζητήσουμε →
              </a>
            </div>
          </div>
          <div className="grid gap-6 sm:grid-cols-2">
            {projects.map((project) => (
              <a
                key={project.title + project.subtitle}
                className="flex h-full flex-col gap-4 rounded-2xl border border-stone-900/10 bg-gradient-to-br from-white to-[#e6f0ff] p-6 transition hover:-translate-y-1 hover:border-[var(--accent)]/50 hover:shadow-lg"
                href={project.href}
                target="_blank"
                rel="noreferrer"
              >
                <div>
                  <p className="text-xs font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">
                    {project.subtitle}
                  </p>
                  <h3 className="mt-2 text-xl font-semibold text-stone-900">
                    {project.title}
                  </h3>
                </div>
                <p className="text-sm text-stone-600">{project.description}</p>
                <span className="mt-auto inline-flex items-center gap-2 text-sm font-semibold text-[var(--accent)]">
                  Δες το project →
                  <span aria-hidden className="block h-2 w-2 rounded-full bg-[var(--accent)]" />
                </span>
              </a>
            ))}
          </div>
        </section>

        <section className="grid gap-10 lg:grid-cols-[1.1fr_0.9fr] lg:items-center">
          <div className="space-y-6">
            <p className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">
              Μεθοδολογία
            </p>
            <h2 className="text-3xl font-semibold sm:text-4xl">
              Ένα συνεργατικό ταξίδι, με διαφάνεια και επίκεντρο τον χρήστη.
            </h2>
            <p className="text-base text-stone-600">
              Συνδυάζουμε στρατηγική, σχεδίαση και ανάπτυξη σε έναν κύκλο ζωής
              που προσαρμόζεται στις ανάγκες σας. Από την πρώτη συνάντηση έως το
              go-live και τα συνεχή iterations, είμαστε η ομάδα που θέλετε στο
              πλευρό σας.
            </p>
          </div>
          <div className="space-y-6 rounded-3xl border border-stone-900/10 bg-white p-8 shadow-sm">
            {milestones.map((milestone, index) => (
              <div key={milestone.title} className="space-y-2 border-l-2 border-[var(--accent)]/40 pl-6">
                <span className="text-xs font-semibold uppercase tracking-[0.3em] text-[var(--accent)]">
                  Βήμα {index + 1}
                </span>
                <h3 className="text-lg font-semibold text-stone-800">
                  {milestone.title}
                </h3>
                <p className="text-sm text-stone-600">{milestone.text}</p>
                {index === 1 ? (
                  <p className="text-xs font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">
                    QA • Monitoring • Green Deploys
                  </p>
                ) : null}
              </div>
            ))}
          </div>
        </section>

        <section className="grid gap-10 rounded-3xl border border-stone-900/10 bg-white p-10 shadow-sm lg:grid-cols-[1.1fr_0.9fr] lg:items-center">
          <div className="space-y-6">
            <p className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">
              Γιατί εμάς
            </p>
            <h2 className="text-3xl font-semibold sm:text-4xl">
              Εμπιστευτείτε μια ομάδα που συνδυάζει δημιουργικότητα και τεχνογνωσία.
            </h2>
            <p className="text-base text-stone-600">
              Κάθε project συνοδεύεται από dedicated project manager, senior
              engineers και designers. Ερευνούμε, σχεδιάζουμε, υλοποιούμε και
              υποστηρίζουμε με πάθος για ποιότητα.
            </p>
          </div>
          <ul className="space-y-4 text-sm text-stone-700">
            <li className="flex items-start gap-3 rounded-2xl border border-stone-900/10 bg-[#eef3ff] p-4">
              <span className="mt-0.5 text-lg text-[var(--accent)]">◆</span>
              <div>
                <p className="font-semibold text-stone-900">Performance-first</p>
                <p className="text-stone-600">
                  Lighthouse 90+ ως προεπιλογή με συνεχές monitoring και
                  βελτιστοποιήσεις SEO.
                </p>
              </div>
            </li>
            <li className="flex items-start gap-3 rounded-2xl border border-stone-900/10 bg-[#eef3ff] p-4">
              <span className="mt-0.5 text-lg text-[var(--accent)]">◆</span>
              <div>
                <p className="font-semibold text-stone-900">Μακροχρόνια σχέση</p>
                <p className="text-stone-600">
                  Συμβόλαια υποστήριξης, SLA και roadmap sessions ανά τρίμηνο.
                </p>
              </div>
            </li>
            <li className="flex items-start gap-3 rounded-2xl border border-stone-900/10 bg-[#eef3ff] p-4">
              <span className="mt-0.5 text-lg text-[var(--accent)]">◆</span>
              <div>
                <p className="font-semibold text-stone-900">Διαφανής κοστολόγηση</p>
                <p className="text-stone-600">
                  Πακέτα και custom προσφορές με πλήρη ανάλυση ώρας και
                  παραδοτέων.
                </p>
              </div>
            </li>
          </ul>
        </section>

        <section className="rounded-3xl border border-stone-900/10 bg-white p-10 shadow-sm">
          <div className="space-y-4">
            <p className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">
              {philosophy.emoji} {philosophy.title}
            </p>
            <h2 className="text-3xl font-semibold text-stone-900 sm:text-4xl">
              Μια μικρή ομάδα με βαθύ χρόνο για τα project της.
            </h2>
            <div className="space-y-3 text-base text-stone-600">
              {philosophy.paragraphs.map((paragraph) => (
                <p key={paragraph}>{paragraph}</p>
              ))}
            </div>
          </div>
        </section>

        <section
          id="contact"
          className="relative overflow-hidden rounded-3xl border border-[var(--accent)]/30 bg-gradient-to-br from-[#e0ebff] via-[#f3f7ff] to-[var(--background)] p-10"
        >
          <div className="absolute inset-0 -z-10 bg-[radial-gradient(circle_at_bottom_right,_rgba(37,99,235,0.18)_0%,_transparent_55%)]" />
          <div className="space-y-6">
            <p className="text-sm font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">
              Επικοινωνία
            </p>
            <h2 className="text-3xl font-semibold sm:text-4xl">
              Ξεκινήστε σήμερα ένα project που θα μιλάει τη γλώσσα της αγοράς σας.
            </h2>
            <p className="max-w-2xl text-base text-stone-700 sm:text-lg">
              Στείλτε μας μια σύντομη περιγραφή των στόχων σας και θα επιστρέψουμε
              με πρόταση, χρονοδιάγραμμα και ενδεικτικό προϋπολογισμό μέσα σε 2
              εργάσιμες ημέρες.
            </p>
            <div className="flex flex-col gap-3 sm:flex-row">
              <a
                className="inline-flex items-center justify-center rounded-full bg-[var(--accent)] px-6 py-3 text-sm font-semibold text-white transition hover:bg-[var(--accent-dark)]"
                href="mailto:hello@outwardly.gr"
              >
                hello@outwardly.gr
              </a>
              <a
                className="inline-flex items-center justify-center rounded-full border border-[var(--accent)] px-6 py-3 text-sm font-semibold text-[var(--accent)] transition hover:bg-[var(--accent)] hover:text-white"
                href="/intro"
              >
                Book intro call →
              </a>
            </div>
          </div>
        </section>
      </main>

      <footer className="border-t border-stone-300 bg-[#e7edff]">
        <div className="mx-auto grid max-w-6xl gap-6 px-6 py-8 text-sm text-stone-600 sm:grid-cols-3 sm:px-8">
          <div className="space-y-1 text-stone-700">
            <p className="font-semibold text-stone-900">© 2025 Outwardly</p>
            <p>ΑΦΜ: [συμπλήρωσέ το]</p>
            <p>ΔΟΥ: Πατρών</p>
          </div>
          <div className="space-y-1">
            <p>ΚΑΔ: 62.01.11.01 – Υπηρεσίες σχεδίασης και ανάπτυξης ιστοσελίδων</p>
            <p>ΚΑΔ: 63.11.11.00 – Υπηρεσίες φιλοξενίας ιστοσελίδων</p>
          </div>
          <div className="space-y-1">
            <a className="block hover:text-[var(--accent)]" href="#services">
              Υπηρεσίες
            </a>
            <a className="block hover:text-[var(--accent)]" href="#contact">
              Επικοινωνία
            </a>
            <a className="block hover:text-[var(--accent)]" href="mailto:hello@outwardly.gr">
              hello@outwardly.gr
            </a>
          </div>
        </div>
      </footer>
    </div>
  );
}
