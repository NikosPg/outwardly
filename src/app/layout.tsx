import type { Metadata } from "next";
import { Geist, Geist_Mono } from "next/font/google";
import Script from "next/script";
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
      <head>
        <Script
          src="https://www.googletagmanager.com/gtag/js?id=G-48EB83B9YY"
          strategy="afterInteractive"
        />
        <Script id="google-analytics" strategy="afterInteractive">
          {`
            window.dataLayer = window.dataLayer || [];
            function gtag(){dataLayer.push(arguments);}
            gtag('js', new Date());
            gtag('config', 'G-48EB83B9YY');
          `}
        </Script>
      </head>
      <body
        className={`${geistSans.variable} ${geistMono.variable} antialiased`}
      >
        {children}
      </body>
    </html>
  );
}
