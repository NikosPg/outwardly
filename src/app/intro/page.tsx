import type { Metadata } from "next";
import { Suspense } from "react";
import { IntroContent } from "./IntroContent";

export const metadata: Metadata = {
  title: "Book an Intro Call • ΕΚΦΑΝΣΙΣ",
  description:
    "Κλείστε ένα σύντομο intro call με την ομάδα της ΕΚΦΑΝΣΙΣ και ξεκινήστε τη συζήτηση για το επόμενο ψηφιακό σας project.",
};

export default function IntroPage({ searchParams }: { searchParams: { lang?: string } }) {
  const defaultLocale = searchParams?.lang === "en" ? "en" : "el";
  return (
    <main className="min-h-screen bg-[var(--background)] text-[var(--foreground)]">
      <Suspense fallback={<div className="min-h-screen" />}>
        <IntroContent defaultLocale={defaultLocale} />
      </Suspense>
    </main>
  );
}
