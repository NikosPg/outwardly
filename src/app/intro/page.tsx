import type { Metadata } from "next";

export const metadata: Metadata = {
  title: "Book an Intro Call • OutWardly",
  description:
    "Κλείστε ένα σύντομο intro call με την ομάδα της OutWardly και ξεκινήστε τη συζήτηση για το επόμενο ψηφιακό σας project.",
};

export default function IntroPage() {
  return (
    <main className="min-h-screen bg-[var(--background)] text-[var(--foreground)]">
      <div className="mx-auto flex max-w-5xl flex-col gap-12 px-6 py-16 sm:px-8 lg:py-24">
        <header className="space-y-4 text-center sm:text-left">
          <p className="text-sm font-semibold uppercase tracking-[0.28em] text-[var(--accent)]">
            Book Intro Call
          </p>
          <h1 className="text-4xl font-semibold sm:text-5xl">
            Ας συστηθούμε και ας θέσουμε τα επόμενα βήματα.
          </h1>
          <p className="text-base text-stone-600 sm:text-lg">
            Επιλέξτε ώρα που σας βολεύει για ένα σύντομο call. Θα μιλήσουμε για τους
            στόχους σας, θα απαντήσουμε σε ερωτήσεις και θα προτείνουμε το πλάνο
            συνεργασίας που ταιριάζει στην ομάδα σας.
          </p>
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
          Αν προτιμάτε email, μπορείτε πάντα να μας γράψετε στο{" "}
          <a className="font-semibold text-[var(--accent)] hover:text-[var(--accent-dark)]" href="mailto:hello@outwardly.gr">
            hello@outwardly.gr
          </a>
          .
        </p>
      </div>
    </main>
  );
}
