"use client";

import { useEffect, useState } from "react";

interface AccessibilitySettings {
  largeText: boolean;
  highContrast: boolean;
  reducedMotion: boolean;
  dyslexiaFont: boolean;
  linkHighlight: boolean;
}

const defaultSettings: AccessibilitySettings = {
  largeText: false,
  highContrast: false,
  reducedMotion: false,
  dyslexiaFont: false,
  linkHighlight: false,
};

export function AccessibilityWidget() {
  const [isOpen, setIsOpen] = useState(false);
  const [settings, setSettings] = useState<AccessibilitySettings>(defaultSettings);

  // Load settings from localStorage on mount
  useEffect(() => {
    const saved = localStorage.getItem("accessibility-settings");
    if (saved) {
      const parsed = JSON.parse(saved) as AccessibilitySettings;
      setSettings(parsed);
      applySettings(parsed);
    }
  }, []);

  const applySettings = (newSettings: AccessibilitySettings) => {
    const root = document.documentElement;

    // Large text
    if (newSettings.largeText) {
      root.style.fontSize = "120%";
    } else {
      root.style.fontSize = "";
    }

    // High contrast
    if (newSettings.highContrast) {
      root.classList.add("high-contrast");
    } else {
      root.classList.remove("high-contrast");
    }

    // Reduced motion
    if (newSettings.reducedMotion) {
      root.classList.add("reduced-motion");
    } else {
      root.classList.remove("reduced-motion");
    }

    // Dyslexia-friendly font
    if (newSettings.dyslexiaFont) {
      root.classList.add("dyslexia-font");
    } else {
      root.classList.remove("dyslexia-font");
    }

    // Highlight links
    if (newSettings.linkHighlight) {
      root.classList.add("highlight-links");
    } else {
      root.classList.remove("highlight-links");
    }
  };

  const updateSetting = (key: keyof AccessibilitySettings) => {
    const newSettings = { ...settings, [key]: !settings[key] };
    setSettings(newSettings);
    localStorage.setItem("accessibility-settings", JSON.stringify(newSettings));
    applySettings(newSettings);
  };

  const resetSettings = () => {
    setSettings(defaultSettings);
    localStorage.removeItem("accessibility-settings");
    applySettings(defaultSettings);
  };

  return (
    <>
      {/* Accessibility Button */}
      <button
        type="button"
        onClick={() => setIsOpen(!isOpen)}
        className="fixed bottom-6 left-6 z-50 flex h-14 w-14 items-center justify-center rounded-full bg-[var(--accent)] text-white shadow-lg transition hover:bg-[var(--accent-dark)] focus:outline-none focus:ring-4 focus:ring-[var(--accent)]/50"
        aria-label="Άνοιγμα ρυθμίσεων προσβασιμότητας"
        aria-expanded={isOpen}
        aria-controls="accessibility-panel"
      >
        {/* Wheelchair/Accessibility Icon */}
        <svg
          xmlns="http://www.w3.org/2000/svg"
          viewBox="0 0 24 24"
          fill="currentColor"
          className="h-7 w-7"
          aria-hidden="true"
        >
          <path d="M12 2a2 2 0 1 1 0 4 2 2 0 0 1 0-4zm8 18h-1l-2.5-5H13v-2h4l1.5 3H20a1 1 0 1 1 0 2zm-9.5-3a3.5 3.5 0 1 1 0-7 3.5 3.5 0 0 1 0 7zm0-2a1.5 1.5 0 1 0 0-3 1.5 1.5 0 0 0 0 3zM9 8v5H7V8a1 1 0 1 1 2 0z"/>
        </svg>
      </button>

      {/* Accessibility Panel */}
      {isOpen && (
        <div
          id="accessibility-panel"
          role="dialog"
          aria-label="Ρυθμίσεις προσβασιμότητας"
          aria-modal="true"
          className="fixed bottom-24 left-6 z-50 w-80 rounded-2xl border border-stone-200 bg-white p-6 shadow-2xl"
        >
          <div className="mb-4 flex items-center justify-between">
            <h2 className="text-lg font-semibold text-stone-900">Προσβασιμότητα</h2>
            <button
              type="button"
              onClick={() => setIsOpen(false)}
              className="rounded-lg p-1 text-stone-500 hover:bg-stone-100 hover:text-stone-700"
              aria-label="Κλείσιμο"
            >
              <svg className="h-5 w-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M6 18L18 6M6 6l12 12" />
              </svg>
            </button>
          </div>

          <div className="space-y-3">
            <ToggleOption
              label="Μεγάλο κείμενο"
              description="Αύξηση μεγέθους γραμματοσειράς"
              checked={settings.largeText}
              onChange={() => updateSetting("largeText")}
              icon={
                <svg className="h-5 w-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M4 6h16M4 12h8m-8 6h16" />
                </svg>
              }
            />

            <ToggleOption
              label="Υψηλή αντίθεση"
              description="Ενισχυμένα χρώματα"
              checked={settings.highContrast}
              onChange={() => updateSetting("highContrast")}
              icon={
                <svg className="h-5 w-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M12 3v1m0 16v1m9-9h-1M4 12H3m15.364 6.364l-.707-.707M6.343 6.343l-.707-.707m12.728 0l-.707.707M6.343 17.657l-.707.707M16 12a4 4 0 11-8 0 4 4 0 018 0z" />
                </svg>
              }
            />

            <ToggleOption
              label="Χωρίς κινήσεις"
              description="Απενεργοποίηση animations"
              checked={settings.reducedMotion}
              onChange={() => updateSetting("reducedMotion")}
              icon={
                <svg className="h-5 w-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 9v6m4-6v6m7-3a9 9 0 11-18 0 9 9 0 0118 0z" />
                </svg>
              }
            />

            <ToggleOption
              label="Γραμματοσειρά δυσλεξίας"
              description="Ευανάγνωστη γραμματοσειρά"
              checked={settings.dyslexiaFont}
              onChange={() => updateSetting("dyslexiaFont")}
              icon={
                <svg className="h-5 w-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 12h6m-6 4h6m2 5H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
                </svg>
              }
            />

            <ToggleOption
              label="Επισήμανση συνδέσμων"
              description="Υπογράμμιση όλων των links"
              checked={settings.linkHighlight}
              onChange={() => updateSetting("linkHighlight")}
              icon={
                <svg className="h-5 w-5" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                  <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M13.828 10.172a4 4 0 00-5.656 0l-4 4a4 4 0 105.656 5.656l1.102-1.101m-.758-4.899a4 4 0 005.656 0l4-4a4 4 0 00-5.656-5.656l-1.1 1.1" />
                </svg>
              }
            />
          </div>

          <button
            type="button"
            onClick={resetSettings}
            className="mt-4 w-full rounded-lg border border-stone-300 px-4 py-2 text-sm font-medium text-stone-700 transition hover:bg-stone-50"
          >
            Επαναφορά ρυθμίσεων
          </button>
        </div>
      )}
    </>
  );
}

interface ToggleOptionProps {
  label: string;
  description: string;
  checked: boolean;
  onChange: () => void;
  icon: React.ReactNode;
}

function ToggleOption({ label, description, checked, onChange, icon }: ToggleOptionProps) {
  return (
    <button
      type="button"
      onClick={onChange}
      className={`flex w-full items-center gap-3 rounded-xl p-3 text-left transition ${
        checked ? "bg-[var(--accent)]/10 ring-2 ring-[var(--accent)]" : "bg-stone-50 hover:bg-stone-100"
      }`}
      role="switch"
      aria-checked={checked}
    >
      <div className={`flex h-10 w-10 items-center justify-center rounded-lg ${checked ? "bg-[var(--accent)] text-white" : "bg-stone-200 text-stone-600"}`}>
        {icon}
      </div>
      <div className="flex-1">
        <p className={`text-sm font-medium ${checked ? "text-[var(--accent)]" : "text-stone-900"}`}>{label}</p>
        <p className="text-xs text-stone-500">{description}</p>
      </div>
      <div className={`h-6 w-11 rounded-full p-1 transition ${checked ? "bg-[var(--accent)]" : "bg-stone-300"}`}>
        <div className={`h-4 w-4 rounded-full bg-white transition-transform ${checked ? "translate-x-5" : "translate-x-0"}`} />
      </div>
    </button>
  );
}
