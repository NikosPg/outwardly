"use client";

import { useState } from "react";

const EMAIL = "hello@outwardly.net";

const translations = {
  el: {
    switchLabel: "EN",
    switchAria: "Switch to English",
    label: "Book Intro Call",
    heading: "Ας συστηθούμε και ας θέσουμε τα επόμενα βήματα.",
    body: "Επιλέξτε ώρα που σας βολεύει για ένα σύντομο call. Θα μιλήσουμε για τους στόχους σας, θα απαντήσουμε σε ερωτήσεις και θα προτείνουμε το πλάνο συνεργασίας που ταιριάζει στην ομάδα σας.",
    emailPrompt: "Αν προτιμάτε email, μπορείτε πάντα να μας γράψετε στο",
  },
  en: {
    switchLabel: "ΕΛ",
    switchAria: "Switch to Greek",
    label: "Book Intro Call",
    heading: "Let’s meet and map the next steps together.",
    body: "Pick a time that suits you for a quick call. We’ll discuss your goals, answer questions, and outline the collaboration plan that fits your team.",
    emailPrompt: "Prefer email? Drop us a line at",
  },
} as const;

type Locale = keyof typeof translations;

export function IntroContent({ defaultLocale }: { defaultLocale: Locale }) {
  const [locale, setLocale] = useState<Locale>(defaultLocale);
  const t = translations[locale];
  const alternateLocale = locale === "el" ? "en" : "el";

  return (
    <div className="relative mx-auto flex max-w-5xl flex-col gap-12 px-6 py-16 sm:px-8 lg:py-24">
      <button
        type="button"
        onClick={() => setLocale(alternateLocale)}
        className="absolute right-6 top-6 inline-flex items-center justify-center rounded-full border border-[var(--accent)]/40 bg-white/80 px-4 py-2 text-xs font-semibold uppercase tracking-[0.28em] text-[var(--accent)] shadow-sm backdrop-blur transition hover:bg-white"
        aria-label={t.switchAria}
      >
        {t.switchLabel}
      </button>
      <header className="space-y-4 text-center sm:text-left">
        <p className="text-sm font-semibold uppercase tracking-[0.28em] text-[var(--accent)]">{t.label}</p>
        <h1 className="text-4xl font-semibold sm:text-5xl">{t.heading}</h1>
        <p className="text-base text-stone-600 sm:text-lg">{t.body}</p>
      </header>

      <section className="overflow-hidden rounded-3xl border border-stone-900/10 bg-white shadow-lg">
        <iframe
          src="https://cal.com/outwardly/intro?embed=1&primaryColor=2563eb"
          title="OutWardly Intro Call Scheduler"
          className="h-[720px] w-full"
          frameBorder="0"
          allow="camera; microphone; fullscreen"
          loading="lazy"
        />
      </section>

      <p className="text-center text-sm text-stone-500">
        {t.emailPrompt}{" "}
        <a className="font-semibold text-[var(--accent)] hover:text-[var(--accent-dark)]" href={`mailto:${EMAIL}`}>
          {EMAIL}
        </a>
        .
      </p>
    </div>
  );
}
