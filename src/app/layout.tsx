import type { Metadata } from "next";
import { Geist, Geist_Mono } from "next/font/google";
import "./globals.css";

const geistSans = Geist({
  variable: "--font-geist-sans",
  subsets: ["latin"],
});

const geistMono = Geist_Mono({
  variable: "--font-geist-mono",
  subsets: ["latin"],
});

export const metadata: Metadata = {
  title: {
    default: "OutWardly – Creative Studio & Custom Software",
    template: "%s • OutWardly",
  },
  description:
    "OutWardly σχεδιάζει ψηφιακές εμπειρίες, websites και custom λογισμικό με focus στην απόδοση και τη στρατηγική.",
  metadataBase: new URL("https://outwardly.net"),
  openGraph: {
    title: "OutWardly – Creative Studio & Custom Software",
    description:
      "Αναλαμβάνουμε end-to-end ανάπτυξη web projects, custom λογισμικό και φιλοξενία για ομάδες που χρειάζονται αξιόπιστη digital παρουσία.",
    url: "https://outwardly.net",
    siteName: "OutWardly",
    locale: "el_GR",
    type: "website",
  },
  twitter: {
    card: "summary_large_image",
    title: "OutWardly – Creative Studio & Custom Software",
    description:
      "Στρατηγική, design, ανάπτυξη και φιλοξενία σε ένα ενιαίο digital studio.",
  },
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body
        className={`${geistSans.variable} ${geistMono.variable} antialiased`}
      >
        {children}
      </body>
    </html>
  );
}
