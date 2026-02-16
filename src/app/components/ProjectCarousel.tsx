"use client";

import Image from "next/image";
import { useState, useEffect, useCallback, useRef } from "react";

interface ProjectSlide {
  title: string;
  subtitle: string;
  description: string;
  cta: string;
  href: string;
  image: string;
}

interface ProjectCarouselProps {
  projects: readonly ProjectSlide[];
  ariaLabels: {
    prev: string;
    next: string;
    slideLabel: string;
  };
}

export function ProjectCarousel({ projects, ariaLabels }: ProjectCarouselProps) {
  const [currentIndex, setCurrentIndex] = useState(0);
  const [isHovered, setIsHovered] = useState(false);
  const [isFocused, setIsFocused] = useState(false);
  const [prefersReducedMotion, setPrefersReducedMotion] = useState(false);
  const pointerStartX = useRef<number | null>(null);
  const pointerStartTime = useRef(0);

  // Detect reduced motion preference (system + manual toggle)
  useEffect(() => {
    const check = () =>
      document.documentElement.classList.contains("reduced-motion") ||
      window.matchMedia("(prefers-reduced-motion: reduce)").matches;

    setPrefersReducedMotion(check());

    const mq = window.matchMedia("(prefers-reduced-motion: reduce)");
    const mqHandler = () => setPrefersReducedMotion(check());
    mq.addEventListener("change", mqHandler);

    const observer = new MutationObserver(() => setPrefersReducedMotion(check()));
    observer.observe(document.documentElement, {
      attributes: true,
      attributeFilter: ["class"],
    });

    return () => {
      mq.removeEventListener("change", mqHandler);
      observer.disconnect();
    };
  }, []);

  const goToNext = useCallback(() => {
    setCurrentIndex((prev) => (prev + 1) % projects.length);
  }, [projects.length]);

  const goToPrev = useCallback(() => {
    setCurrentIndex((prev) => (prev - 1 + projects.length) % projects.length);
  }, [projects.length]);

  // Auto-advance every 5 seconds
  useEffect(() => {
    if (isHovered || isFocused || prefersReducedMotion) return;
    const timer = setInterval(goToNext, 5000);
    return () => clearInterval(timer);
  }, [currentIndex, isHovered, isFocused, prefersReducedMotion, goToNext]);

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === "ArrowLeft") {
      e.preventDefault();
      goToPrev();
    } else if (e.key === "ArrowRight") {
      e.preventDefault();
      goToNext();
    }
  };

  const handlePointerDown = (e: React.PointerEvent) => {
    pointerStartX.current = e.clientX;
    pointerStartTime.current = Date.now();
    (e.target as HTMLElement).setPointerCapture?.(e.pointerId);
  };

  const handlePointerUp = (e: React.PointerEvent) => {
    if (pointerStartX.current === null) return;
    const deltaX = e.clientX - pointerStartX.current;
    const elapsed = Date.now() - pointerStartTime.current;
    const velocity = Math.abs(deltaX) / elapsed;

    if (Math.abs(deltaX) > 50 || velocity > 0.3) {
      if (deltaX < 0) goToNext();
      else goToPrev();
    }
    pointerStartX.current = null;
  };

  const slideLabel = ariaLabels.slideLabel
    .replace("{current}", String(currentIndex + 1))
    .replace("{total}", String(projects.length));

  return (
    <div
      role="region"
      aria-roledescription="carousel"
      aria-label={slideLabel}
      onMouseEnter={() => setIsHovered(true)}
      onMouseLeave={() => setIsHovered(false)}
      onFocus={() => setIsFocused(true)}
      onBlur={(e) => {
        if (!e.currentTarget.contains(e.relatedTarget as Node)) {
          setIsFocused(false);
        }
      }}
      onKeyDown={handleKeyDown}
      tabIndex={0}
      className="relative outline-none focus-visible:outline-3 focus-visible:outline-[var(--accent)] focus-visible:rounded-2xl"
    >
      {/* Screen reader live region */}
      <div aria-live="polite" aria-atomic="true" className="sr-only">
        {slideLabel}: {projects[currentIndex].title}
      </div>

      {/* Slides viewport */}
      <div
        className="overflow-hidden rounded-2xl"
        style={{ touchAction: "pan-y" }}
        onPointerDown={handlePointerDown}
        onPointerUp={handlePointerUp}
      >
        <div
          className="flex transition-transform duration-500 ease-in-out"
          style={{
            transform: `translateX(-${currentIndex * 100}%)`,
            transitionDuration: prefersReducedMotion ? "0ms" : "500ms",
          }}
        >
          {projects.map((project, index) => (
            <div
              key={project.title}
              role="group"
              aria-roledescription="slide"
              aria-label={`${index + 1} / ${projects.length}: ${project.title}`}
              aria-hidden={index !== currentIndex}
              className="w-full flex-shrink-0"
            >
              <div className="grid gap-6 rounded-2xl border border-stone-900/10 bg-white p-5 sm:p-6 lg:grid-cols-2 lg:gap-10 lg:p-8">
                {/* Screenshot */}
                <div className="overflow-hidden rounded-xl border border-stone-900/10 shadow-md">
                  <Image
                    src={project.image}
                    alt={`Screenshot: ${project.title}`}
                    width={1280}
                    height={800}
                    sizes="(max-width: 1024px) 100vw, 50vw"
                    className="h-auto w-full object-cover"
                    priority={index === 0}
                    loading={index === 0 ? "eager" : "lazy"}
                  />
                </div>

                {/* Text content */}
                <div className="flex flex-col justify-center gap-3 sm:gap-4">
                  <p className="text-xs font-semibold uppercase tracking-[0.2em] text-[var(--accent)]">
                    {project.subtitle}
                  </p>
                  <h3 className="text-2xl font-semibold text-stone-900 sm:text-3xl">
                    {project.title}
                  </h3>
                  <p className="text-sm leading-relaxed text-stone-600 sm:text-base">
                    {project.description}
                  </p>
                  <a
                    href={project.href}
                    target="_blank"
                    rel="noreferrer"
                    tabIndex={index === currentIndex ? 0 : -1}
                    className="mt-2 inline-flex items-center gap-2 text-sm font-semibold text-[var(--accent)] transition hover:text-[var(--accent-dark)]"
                  >
                    {project.cta}
                    <svg
                      className="h-4 w-4"
                      fill="none"
                      viewBox="0 0 24 24"
                      stroke="currentColor"
                      strokeWidth={2}
                      aria-hidden="true"
                    >
                      <path
                        strokeLinecap="round"
                        strokeLinejoin="round"
                        d="M14 5l7 7m0 0l-7 7m7-7H3"
                      />
                    </svg>
                  </a>
                </div>
              </div>
            </div>
          ))}
        </div>
      </div>

      {/* Left arrow */}
      <button
        type="button"
        onClick={goToPrev}
        aria-label={ariaLabels.prev}
        className="absolute left-2 top-1/2 z-10 -translate-y-1/2 flex h-10 w-10 items-center justify-center rounded-full border border-stone-900/10 bg-white/90 text-stone-700 shadow-sm backdrop-blur transition hover:bg-white hover:text-[var(--accent)] focus-visible:outline focus-visible:outline-3 focus-visible:outline-[var(--accent)] sm:left-3 sm:h-12 sm:w-12"
      >
        <svg
          className="h-5 w-5"
          fill="none"
          viewBox="0 0 24 24"
          stroke="currentColor"
          strokeWidth={2}
          aria-hidden="true"
        >
          <path strokeLinecap="round" strokeLinejoin="round" d="M15 19l-7-7 7-7" />
        </svg>
      </button>

      {/* Right arrow */}
      <button
        type="button"
        onClick={goToNext}
        aria-label={ariaLabels.next}
        className="absolute right-2 top-1/2 z-10 -translate-y-1/2 flex h-10 w-10 items-center justify-center rounded-full border border-stone-900/10 bg-white/90 text-stone-700 shadow-sm backdrop-blur transition hover:bg-white hover:text-[var(--accent)] focus-visible:outline focus-visible:outline-3 focus-visible:outline-[var(--accent)] sm:right-3 sm:h-12 sm:w-12"
      >
        <svg
          className="h-5 w-5"
          fill="none"
          viewBox="0 0 24 24"
          stroke="currentColor"
          strokeWidth={2}
          aria-hidden="true"
        >
          <path strokeLinecap="round" strokeLinejoin="round" d="M9 5l7 7-7 7" />
        </svg>
      </button>

      {/* Dot indicators */}
      <div className="mt-6 flex justify-center gap-2" role="tablist">
        {projects.map((project, index) => (
          <button
            key={project.title}
            role="tab"
            aria-selected={index === currentIndex}
            aria-label={project.title}
            onClick={() => setCurrentIndex(index)}
            className={`h-2.5 rounded-full transition-all duration-300 ${
              index === currentIndex
                ? "w-8 bg-[var(--accent)]"
                : "w-2.5 bg-stone-300 hover:bg-stone-400"
            }`}
          />
        ))}
      </div>
    </div>
  );
}
