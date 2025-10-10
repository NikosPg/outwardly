import type { Metadata } from "next";
import { IntroContent } from "./IntroContent";

export const metadata: Metadata = {
  title: "Book an Intro Call • OutWardly",
  description:
    "Κλείστε ένα σύντομο intro call με την ομάδα της OutWardly και ξεκινήστε τη συζήτηση για το επόμενο ψηφιακό σας project.",
};

export default function IntroPage({ searchParams }: { searchParams: { lang?: string } }) {
  const defaultLocale = searchParams?.lang === "en" ? "en" : "el";
  return (
    <main className="min-h-screen bg-[var(--background)] text-[var(--foreground)]">
      <IntroContent defaultLocale={defaultLocale} />
    </main>
  );
}
