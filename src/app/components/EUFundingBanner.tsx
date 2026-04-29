export function EUFundingBanner() {
  return (
    <section
      aria-label="Χρηματοδότηση από την Ευρωπαϊκή Ένωση"
      className="border-t border-stone-200 bg-[#f4f6ff] py-8"
    >
      <div className="mx-auto max-w-6xl px-6 sm:px-8">
        <a
          href="https://greece20.gov.gr"
          target="_blank"
          rel="noopener noreferrer"
          className="flex flex-col items-center gap-6 sm:flex-row sm:items-center sm:gap-8"
          aria-label="Ελλάδα 2.0 – Εθνικό Σχέδιο Ανάκαμψης και Ανθεκτικότητας (ανοίγει σε νέα καρτέλα)"
        >
          <div className="flex-shrink-0">
            <img
              src="/greece20-logo.png"
              alt="Ελλάδα 2.0 – Εθνικό Σχέδιο Ανάκαμψης και Ανθεκτικότητας, χρηματοδοτούμενο από την Ευρωπαϊκή Ένωση NextGenerationEU"
              className="h-16 w-auto"
            />
          </div>

          <div className="hidden h-14 w-px bg-stone-300 sm:block" aria-hidden="true" />

          <div className="flex-1 space-y-1 text-center sm:text-left">
            <p className="text-sm font-semibold text-stone-800">
              Με τη χρηματοδότηση της Ευρωπαϊκής Ένωσης – NextGenerationEU
            </p>
            <p className="text-xs text-stone-600">
              Το παρόν έργο υλοποιείται στο πλαίσιο του Εθνικού Σχεδίου Ανάκαμψης και
              Ανθεκτικότητας «Ελλάδα 2.0», με τη χρηματοδότηση της Ευρωπαϊκής Ένωσης –
              NextGenerationEU.
            </p>
          </div>
        </a>
      </div>
    </section>
  );
}
