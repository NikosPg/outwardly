"use client";

import Link from "next/link";
import { useCallback, useEffect, useMemo, useState } from "react";
import { usePathname, useRouter, useSearchParams } from "next/navigation";

const EMAIL = "ekfansis@gmail.com";

const translations = {
  el: {
    label: "Book Intro Call",
    heading: "Ας συστηθούμε και ας θέσουμε τα επόμενα βήματα.",
    body: "Επιλέξτε ώρα που σας βολεύει για ένα σύντομο call. Θα μιλήσουμε για τους στόχους σας, θα απαντήσουμε σε ερωτήσεις και θα προτείνουμε το πλάνο συνεργασίας που ταιριάζει στην ομάδα σας.",
    emailPrompt: "Αν προτιμάτε email, μπορείτε πάντα να μας γράψετε στο",
    backHome: "← Πίσω στην αρχική",
  },
  en: {
    label: "Book Intro Call",
    heading: "Let’s meet and map the next steps together.",
    body: "Pick a time that suits you for a quick call. We’ll discuss your goals, answer questions, and outline the collaboration plan that fits your team.",
    emailPrompt: "Prefer email? Drop us a line at",
    backHome: "← Back to homepage",
  },
} as const;

type Locale = keyof typeof translations;

const LOCALE_SWITCH: Record<Locale, { label: string; aria: string; next: Locale }> = {
  el: { label: "EN", aria: "Switch to English", next: "en" },
  en: { label: "ΕΛ", aria: "Switch to Greek", next: "el" },
};

export function IntroContent({ defaultLocale }: { defaultLocale: Locale }) {
  const searchParams = useSearchParams();
  const router = useRouter();
  const pathname = usePathname();

  const [locale, setLocaleState] = useState<Locale>(defaultLocale);

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

  const homeHref = useMemo(() => (locale === "en" ? "/?lang=en" : "/"), [locale]);

  return (
    <div className="relative mx-auto flex max-w-5xl flex-col gap-12 px-6 py-16 sm:px-8 lg:py-24">
      <div className="absolute left-6 top-6 flex items-center gap-3">
        <Link
          href={homeHref}
          className="inline-flex items-center rounded-full border border-[var(--accent)]/40 bg-white/80 px-4 py-2 text-xs font-semibold uppercase tracking-[0.28em] text-[var(--accent)] shadow-sm backdrop-blur transition hover:bg-white"
        >
          {t.backHome}
        </Link>
        <button
          type="button"
          onClick={() => handleLocaleChange(LOCALE_SWITCH[locale].next)}
          className="inline-flex items-center rounded-full border border-[var(--accent)]/40 bg-white/80 px-4 py-2 text-xs font-semibold uppercase tracking-[0.28em] text-[var(--accent)] shadow-sm backdrop-blur transition hover:bg-white"
          aria-label={LOCALE_SWITCH[locale].aria}
        >
          {LOCALE_SWITCH[locale].label}
        </button>
      </div>

      <header className="mt-16 space-y-4 text-center sm:mt-0 sm:text-left">
        <p className="text-sm font-semibold uppercase tracking-[0.28em] text-[var(--accent)]">{t.label}</p>
        <h1 className="text-4xl font-semibold sm:text-5xl">{t.heading}</h1>
        <p className="text-base text-stone-600 sm:text-lg">{t.body}</p>
      </header>

      <section className="flex flex-col items-center justify-center rounded-3xl border border-stone-900/10 bg-white p-12 shadow-lg">
        <p className="mb-6 text-center text-lg text-stone-600">
          {t.emailPrompt}
        </p>
        <a
          className="inline-flex items-center justify-center rounded-full bg-[var(--accent)] px-8 py-4 text-base font-semibold text-white transition hover:bg-[var(--accent-dark)]"
          href={`mailto:${EMAIL}`}
        >
          {EMAIL}
        </a>
      </section>
    </div>
  );
}
